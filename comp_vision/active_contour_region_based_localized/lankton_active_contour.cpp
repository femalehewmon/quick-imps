#include "characterization.hpp"
#include "chan_vese_active_contour.hpp"

// 2008 localized active contour implementation
//
// Lankton, Shawn, and Allen Tannenbaum.
// "Localizing region-based active contours."
// IEEE transactions on image processing 17.11 (2008): 2029-2039.

// assign circle pattern
void initializePhiCircle(Mat& phi, int width, int height) {
    int w = floor(width/2) - 50;
    int h = floor(height/2) - 50;
    for ( int i = 0; i < phi.rows; ++i ) {
        for ( int j = 0; j < phi.cols; ++j ) {
            phi.at<double>(i, j) = -sqrt(pow(j-w-50, 2) + pow(i-h-50, 2)) + min(w, h);
        }
    }
}

double areaInsideContour(Mat phi) {
    return (double)foreground<double>(phi).size();
}



double diracDelta (double phi_n, double dt) {
    // d/dt of H(phi) = 1/2*(1 + (2/PI)*arctan(phi))
    return dt / (M_PI*(pow(dt, 2) + pow(phi_n, 2)));
}

//
// Computes a single iteration of the Lankton active contour algorithm.
// Returns true if the contour converged, false otherwise.
//
bool activeContour(Mat src, Mat& phi) {

    double lambda = 1.0;        // contour arc length weight
    double radius = 10;         // radius of localized region to consider

    double dt = 0.25;                   // artificial time step t >= 0
    double threshold_tolerance = 0.01;  // stopping condition

    double area = areaInsideContour(phi);
    double c1 = regionAverage(src, phi, INSIDE);  // pixel avg inside contour
    double c2 = regionAverage(src, phi, OUTSIDE); // pixel avg outside contour
    cout << "c1: " << c1 << " c2: " << c2 << endl;

    int u0 = 0;
    double phi_n, phi_n_plus_1 = 0.0;
    double delta = 0.0;
    double total_energy;
    Point p;

    double phi_diff = 0.0;

    int B = 0;
	int px_i, px_j;
	int py_i, py_j;
    // TODO: fix edge cases, currently not processing borders of image
    for (x_i = 1; x_i < phi.rows-1; ++x_i) {
        for (int x_j = 1; x_j < phi.cols-1; ++x_j) {
            pointX = Point(x_j, x_i);                // current point
            delta = diracDelta(phi_n, dt);           // dirac delta

            if ( delta > 0 ) {
                local_px_i0 = max(0, x_i - radius);
                local_px_in = min(phi.rows, x_i + radius);
                local_px_j0 = max(0, x_j - radius);
                local_px_jn = min(phi.cols, x_j + radius);
                // calculate energy within local region surrounding contour point
                // traverse square only to minimize process time
                for (y_i = local_px_i0; y_i < local_px_in; ++y_i) {
                    for (y_j = local_px_j0; y_j < local_px_jn; ++y_j) {
                        pointY = Point(y_j, y_i);           // prospective point
                        B = withinLocalRegion(pointX, pointY, radius);

                    }
                }
            }
            u0 = src.at<uchar>(i, j);           // current point intensity
            phi_n = phi.at<double>(i, j);        // current phi value

            // Calculate arc length energy
            // --------------------------------------
            double phix, phiy, divr, divl, divu, divd = 0.0;
            double pcurr = phi.at<double>(p.y, p.x);
            phix = phi.at<double>(p.y, p.x+1) - pcurr;
            phiy = (phi.at<double>(p.y+1, p.x) - phi.at<double>(p.y-1, p.x)) / 2;
            divr = (1/sqrt(h + pow(phix, 2) + pow(phiy, 2)));
            // ------------
            phix = pcurr - phi.at<double>(p.y, p.x-1);
            divl = (1/sqrt(h + pow(phix, 2) + pow(phiy, 2)));
            // ------------
            phix = (phi.at<double>(p.y, p.x+1) - phi.at<double>(p.y, p.x-1)) / 2;
            phiy = phi.at<double>(p.y+1, p.x) - pcurr;
            divd = (1/sqrt(h + pow(phix, 2) + pow(phiy, 2)));
            // ------------
            phiy = pcurr - phi.at<double>(p.y-1, p.x);
            divu = (1/sqrt(h + pow(phix, 2) + pow(phiy, 2)));
            // --------------------------------------

            double dist1 = pow(u0 - c1, 2);
            double dist2 = pow(u0 - c2, 2);

            double div = (1 + dt*delta*mu*(divr + divl + divd + divu));
            total_energy = (pcurr + dt*delta*(
                        mu*(phi.at<double>(p.y,p.x+1)*divr +
                            phi.at<double>(p.y,p.x-1)*divl +
                            phi.at<double>(p.y+1,p.x)*divd +
                            phi.at<double>(p.y-1,p.x)*divu)
                        - v*area - lambda1*dist1 + lambda2*dist2)) / div;

            phi_diff += abs(pcurr - total_energy);
            phi.at<double>(i, j) = total_energy;
        }
    }

    // normalize total difference in old vs new phi (describes the contour)
    phi_diff = phi_diff / (phi.rows * phi.cols);
    cout << "difference: " << phi_diff << endl;

    // TODO: optionally reinitialize phi

    return (phi_diff < threshold_tolerance);
}

void getContourMaskFromPhi(Mat phi, Mat& contourMask) {
    contourMask = Mat::zeros(phi.size(), CV_8UC1);
    vector<Point> contour = foreground<double>(phi);
    for (auto point: contour) {
        contourMask.at<uchar>(point.y, point.x) = 255;
    }
}

