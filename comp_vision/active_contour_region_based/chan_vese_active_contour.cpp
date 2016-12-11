#include "characterization.hpp"
#include "chan_vese_active_contour.hpp"

// 2001 Chan-Vese active contour implementation
//
// Chan, Tony F., and Luminita A. Vese.
//      "Active contours without edges."
//      IEEE Transactions on image processing 10.2 (2001): 266-277.
//

// assign checkerboard pattern from IPOL
void initializePhiCheckerboard(Mat& phi) {
    for ( int i = 0; i < phi.rows; ++i ) {
        for ( int j = 0; j < phi.cols; ++j ) {
            phi.at<double>(i, j) = sin(M_PI/5*i)*sin(M_PI/5*j);
        }
    }
}

// assign circle pattern
void initializePhiCircle(Mat& phi, int width, int height) {
    int w = floor(width/2);
    int h = floor(height/2);
    for ( int i = 0; i < phi.rows; ++i ) {
        for ( int j = 0; j < phi.cols; ++j ) {
            phi.at<double>(i, j) = -sqrt(pow(j-w, 2) + pow(i-h, 2)) + min(w, h);
        }
    }
}

int areaInsideContour(Mat phi) {
    int area = 0;
    for ( int row = 0; row < phi.rows; ++row ) {
        for ( int col = 0; col < phi.cols; ++col ) {
            if (phi.at<double>(row, col) >= 0) {
                ++area;
            }
        }
    }
    return area;
}

double regionAverage(Mat src, Mat phi, RegionType region) {
    assert (src.size() == phi.size());
    double total_intensity = 0.0;
    int area = 0;
    for (int row = 0; row < phi.rows; ++row) {
        for ( int col = 0; col < phi.cols; ++col ) {
            if ((region == INSIDE && phi.at<double>(row, col) >= 0) ||
                (region == OUTSIDE && phi.at<double>(row, col) < 0)) {
                ++area;
                total_intensity += src.at<uchar>(row, col);
            }
        }
    }
    return total_intensity / area;
}

double diracDelta (double phi_n, double dt) {
    // d/dt of H(phi) = 1/2*(1 + (2/PI)*arctan(phi))
    return dt / M_PI*(1 + pow(phi_n, 2));
    //return dt / M_PI*(1 + pow(phi_n, 2));
    //return (abs(phi_n - 0.0) < 0.1) ? 1.0 : 0.0;
}

#define DIVIDE_EPS       (1e-16)

double curveEnergy (Mat phi, Point p, double h, double mu) {
    double pcurr = phi.at<double>(p.y, p.x);
    double p_dx1 = pcurr - phi.at<double>(p.y, p.x-1);
    double p_dx2 = phi.at<double>(p.y, p.x+1) - pcurr;
    double p_dx3 = phi.at<double>(p.y, p.x+1) - phi.at<double>(p.y, p.x-1);
    double p_dy1 = pcurr - phi.at<double>(p.y-1, p.x);
    double p_dy2 = phi.at<double>(p.y+1, p.x) - pcurr;
    double p_dy3 = phi.at<double>(p.y+1, p.x) - phi.at<double>(p.y-1, p.x);
    /*
    cout << pcurr << " pcurr" << endl;
    cout << phi.at<double>(p.y, p.x-1) << " pprevx" << endl;
    cout << p_dx1 << " pdx1 " << p_dx2 << " pdx2 " << p_dx3 << " pdx3" << endl;
    cout << p_dy1 << " pdy1 " << p_dy2 << " pdy2 " << p_dy3 << " pdy3" << endl;
    */

    double h1 = pow(h, 2);
    double h2 = pow(2*h, 2);

    double e1 = (mu/h1)*p_dx1
                * sqrt((pow(p_dx2, 2) / h1) + (pow(p_dy3, 2) / h2));
    double e2 = (mu/h1)*p_dy1
                * sqrt((pow(p_dx3, 2) / h2) + (pow(p_dy2, 2) / h1));
    return e1 + e2;
}

//
// Computes a single iteration of the Chan Vese active contour algorithm.
// Returns true if the contour converged, false otherwise.
//
bool activeContour(Mat src, Mat& phi) {

    double mu = 1.0;         // length weight
    double v = 0.0;         // area weight, usually 0 (?)
    double lambda1 = 1.0;   // c1 weight (average intensity inside contour)
    double lambda2 = 1.0;   // c2 weight (average intensity outside contour)

    double h  = 1.0;                    // space step, where (i*h, j*h) are
                                        // the pixel locations
    double dt = 0.5;                    // artificial time step t >= 0
    double threshold_tolerance = 0.01;  // stopping condition

    Mat phi_updated = phi.clone();     // buffer to hold updated phi

    double area = areaInsideContour(phi);
    double c1 = regionAverage(src, phi, INSIDE);  // pixel avg inside contour
    double c2 = regionAverage(src, phi, OUTSIDE); // pixel avg outside contour
    cout << "c1: " << c1 << " c2: " << c2 << endl;

    int u0 = 0;
    double phi_n, phi_n_plus_1 = 0.0;
    double delta = 0.0;
    double total_energy;
    double nu = 0.0;
    Point p;

    double phi_diff = 0.0;
    // TODO: fix edge cases, currently not processing borders of image
    for (int i = 1; i < phi.rows-1; ++i) {
        for (int j = 1; j < phi.cols-1; ++j) {
            u0 = src.at<uchar>(i, j);           // current point intensity
            phi_n = phi.at<double>(i, j);        // current phi value
            delta = diracDelta(phi_n, dt);      // dirac delta
            p = Point(j, i);

            // --------------------------------------
            double phix, phiy, divr, divl, divu, divd = 0.0;
            double pcurr = phi.at<double>(p.y, p.x);
            phix = phi.at<double>(p.y, p.x+1) - pcurr;
            phiy = (phi.at<double>(p.y+1, p.x) - phi.at<double>(p.y-1, p.x)) / 2;
            divr = (1/sqrt(DIVIDE_EPS + pow(phix, 2) + pow(phiy, 2)));
            // ------------
            phix = pcurr - phi.at<double>(p.y, p.x-1);
            divl = (1/sqrt(DIVIDE_EPS + pow(phix, 2) + pow(phiy, 2)));
            // ------------
            phix = (phi.at<double>(p.y, p.x+1) - phi.at<double>(p.y, p.x-1)) / 2;
            phiy = phi.at<double>(p.y+1, p.x) - pcurr;
            divd = (1/sqrt(DIVIDE_EPS + pow(phix, 2) + pow(phiy, 2)));
            // ------------
            phiy = pcurr - phi.at<double>(p.y-1, p.x);
            divu = (1/sqrt(DIVIDE_EPS + pow(phix, 2) + pow(phiy, 2)));
            // --------------------------------------

            double dist1 = abs(u0 - c1);
            double dist2 = abs(u0 - c2);
            //double dist1 = pow(u0 - c1, 2);
            //double dist2 = pow(u0 - c2, 2);

            double div = (1 + delta*mu*(divr + divl + divd + divu));
            total_energy = (div == 0) ? 0 : (pcurr + delta*(mu*(
                            phi.at<double>(p.y,p.x+1)*divr +
                            phi.at<double>(p.y,p.x-1)*divl +
                            phi.at<double>(p.y+1,p.x)*divd +
                            phi.at<double>(p.y-1,p.x)*divu) -
                        nu - lambda1*dist1 + lambda2*dist2)) / div;

            if ( total_energy != total_energy ) {
                cout << "delta " << delta << endl;
                cout << "pcurr " << pcurr << endl;
                cout << "pleft " << phi.at<double>(p.y, p.x-1) << endl;
                cout << "pright " << phi.at<double>(p.y, p.x+1) << endl;
                cout << "pup " << phi.at<double>(p.y-1, p.x) << endl;
                cout << "pdown " << phi.at<double>(p.y+1, p.x) << endl;
                cout << "divr " << divr << endl;
                cout << "divl " << divl << endl;
                cout << "divu " << divu << endl;
                cout << "divd " << divd << endl;
                cout << "dist1 " << dist1 << endl;
                cout << "dist2 " << dist2 << endl;
            }

            total_energy = (total_energy != total_energy) ? 0 : total_energy;

            //phi_diff += pow(pcurr - total_energy, 2);
            phi_diff += abs(pcurr - total_energy);
            phi.at<double>(i, j) = total_energy;

            /*
            // calculate total energy
            total_energy = (curveEnergy(phi, Point(j, i), h, mu) + v*area
                            + lambda1*pow(u0 - c1, 2)
                            + lambda2*pow(u0 - c2, 2)) * delta;
            phi_n_plus_1 = total_energy*dt + phi_n;

            // get difference
            phi_diff += abs(phi.at<double>(i, j) - phi_n_plus_1);

            // update phi buffer
            phi.at<double>(i, j) = phi_n_plus_1;
            */
        }
        if ( total_energy != total_energy ) {
            cout << "\nnan found...\n" << endl;
            break;
        }
    }
    phi_diff = phi_diff / (phi.rows * phi.cols);

    // calculate total difference in old vs new phi (describes the contour)
    //double phi_diff = sum(abs(phi_updated - phi))[0];
    cout << "difference: " << phi_diff << endl;

    // optionally reinitialize phi
    // phi = phi_updated;

    return (phi_diff < threshold_tolerance);
}

void getContourFromPhi(Mat phi, vector<Point>& contour) {
    contour = foreground<double>(phi, 0);
}

