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
            phi.at<float>(i, j) = sin(M_PI/5*i)*sin(M_PI/5*j);
        }
    }
}

// assign circle pattern
void initializePhiCircle(Mat& phi, int width, int height) {
    int w = floor(width/2);
    int h = floor(height/2);
    for ( int i = 0; i < phi.rows; ++i ) {
        for ( int j = 0; j < phi.cols; ++j ) {
            phi.at<float>(i, j) = -sqrt(pow(j-w, 2) + pow(i-h, 2)) + min(w, h);
        }
    }
}

int areaInsideContour(Mat phi) {
    int area = 0;
    for ( int row = 0; row < phi.rows; ++row ) {
        for ( int col = 0; col < phi.cols; ++col ) {
            if (phi.at<float>(row, col) >= 0) {
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
            if ((region == INSIDE && phi.at<float>(row, col) >= 0) ||
                (region == OUTSIDE && phi.at<float>(row, col) < 0)) {
                ++area;
                total_intensity += src.at<uchar>(row, col);
            }
        }
    }
    return total_intensity / area;
}

double diracDelta (double phi_n, double dt) {
    // d/dt of H(phi) = 1/2*(1 + (2/PI)*arctan(phi))
    //return dt / M_PI*(1 + pow(phi_n, 2));
    return (abs(phi_n - 0.0) < 0.1) ? 1.0 : 0.0;
}

double curveEnergyIpol (Mat phi, Point p, double h, double mu) {
    double Aij = mu / sqrt(1);
    double Bij = mu / sqrt(1);
}

double curveEnergy (Mat phi, Point p, double h, double mu) {
    float pcurr = phi.at<float>(p.y, p.x);
    float p_dx1 = pcurr - phi.at<float>(p.y, p.x-1);
    float p_dx2 = phi.at<float>(p.y, p.x+1) - pcurr;
    float p_dx3 = phi.at<float>(p.y, p.x+1) - phi.at<float>(p.y, p.x-1);
    float p_dy1 = pcurr - phi.at<float>(p.y-1, p.x);
    float p_dy2 = phi.at<float>(p.y+1, p.x) - pcurr;
    float p_dy3 = phi.at<float>(p.y+1, p.x) - phi.at<float>(p.y-1, p.x);
    /*
    cout << pcurr << " pcurr" << endl;
    cout << phi.at<float>(p.y, p.x-1) << " pprevx" << endl;
    cout << p_dx1 << " pdx1 " << p_dx2 << " pdx2 " << p_dx3 << " pdx3" << endl;
    cout << p_dy1 << " pdy1 " << p_dy2 << " pdy2 " << p_dy3 << " pdy3" << endl;
    */

    double h1 = pow(h, 2);
    double h2 = pow(2*h, 2);

    float e1 = (mu/h1)*p_dx1
                * sqrt((pow(p_dx2, 2) / h1) + (pow(p_dy3, 2) / h2));
    float e2 = (mu/h1)*p_dy1
                * sqrt((pow(p_dx3, 2) / h2) + (pow(p_dy2, 2) / h1));
    return e1 + e2;
}

//
// Computes a single iteration of the Chan Vese active contour algorithm.
// Returns true if the contour converged, false otherwise.
//
bool activeContour(Mat src, Mat& phi) {

    double u = 1.0;         // length weight
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
    float phi_n, phi_n_plus_1 = 0.0;
    double delta = 0.0;
    double total_energy;

    // TODO: fix edge cases, currently not processing borders of image
    for (int i = 1; i < phi.rows-1; ++i) {
        for (int j = 1; j < phi.cols-1; ++j) {
            u0 = src.at<uchar>(i, j);           // current point intensity
            phi_n = phi.at<float>(i, j);        // current phi value
            delta = diracDelta(phi_n, dt);      // dirac delta

            // calculate total energy
            total_energy = (curveEnergy(phi, Point(j, i), h, u) + v*area
                            + lambda1*pow(u0 - c1, 2)
                            + lambda2*pow(u0 - c2, 2)) * delta;
            phi_n_plus_1 = total_energy*dt + phi_n;
            if ( phi_n_plus_1 == 0 || phi_n == 0 )
                cout << phi_n_plus_1 << " old " << phi_n << endl;
            // update phi buffer
            phi.at<float>(i, j) = phi_n_plus_1;
        }
    }

    // calculate total difference in old vs new phi (describes the contour)
    double phi_diff = sum(abs(phi_updated - phi))[0];
    cout << "difference: " << phi_diff << endl;

    // optionally reinitialize phi (by not doing this)
    //phi = phi_updated.clone();

    return (phi_diff < threshold_tolerance);
}

void getContourFromPhi(Mat phi, vector<Point>& contour) {
    contour = foreground<float>(phi, 0);
}

