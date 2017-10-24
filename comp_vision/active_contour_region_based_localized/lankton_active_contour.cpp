#include "utils.hpp"
#include "characterization.hpp"
#include "lankton_active_contour.hpp"

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

double localRegionAverage(
        Mat src, Mat phi, vector<Point> local_points,RegionType region) {
    assert (local_points.size() > 0);
    double total_intensity = 0.0;
    int area = 0;
    for (auto point: local_points) {
        if ((region == INSIDE && phi.at<double>(point) >= 0) ||
            (region == OUTSIDE && phi.at<double>(point) < 0)) {
            ++area;
            total_intensity += src.at<uchar>(point);
        }
    }
    return area > 0 ? total_intensity / area : 0.0;
}


vector<Point> localSurroundingRegion(Mat image, Point p, int radius) {
    vector<Point> localPoints;
    Point point;
    // scan square quadrant region, save points within circle radius bounds
    for (int i = p.y - radius; i <= p.y; i++) {
        for (int j = p.x - radius; j <= p.x; ++j) {
            point = Point(j, i);           // prospective point
            if (inBounds(image, point) &&
                    distanceBetweenPoints(p, point) <= radius) {
                localPoints.push_back(point);
            }
        }
    }
    return localPoints;
}

/*
vector<Point> getBorderPoints(Mat phi) {
    vector<Point> border_points;
    for(int i = 0; i < phi.rows; ++i) {
        for(int j = 0; j < phi.cols; ++j) {
            if (phi.at<double>(i, j) >= -1.2 &&
                    phi.at<double>(i, j) <= 1.2) {
                border_points.push_back(Point(j, i));
            }
        }
    }
    Mat mask = Mat::zeros(phi.size(), phi.type());
    for (auto p : border_points) {
        circle(mask, p, 2, Scalar(255,0,0), 1);
    }
    imshow("tmp", mask);
    waitKey(0);
    return border_points;
}
*/
vector<Point> getBorderPoints(Mat phi) {
    Mat mask;
    getContourMaskFromPhi(phi, mask);

	Mat edges;
    vector<vector<Point> > contours;
    vector<Vec4i> hierarchy;
    Canny(mask, edges, 1, 100, 3);
    findContours(edges, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE);

    vector<Point> border_points;
    for( int i = 0; i < contours.size(); i++ ) {
        border_points.insert(
                border_points.end(), contours[i].begin(), contours[i].end());
    }

    /*
    mask = Mat::zeros(mask.size(), mask.type());
    for (auto p : border_points) {
        circle(mask, p, 2, Scalar(255,0,0), -1);
    }
    imshow("tmp", mask);
    waitKey(0);
    */

    return border_points;
}

double getCurvature(Mat phi, Point p, double eps) {
    Point above = Point(p.x, p.y+1);
    Point below = Point(p.x, p.y-1);
    Point left  = Point(p.x-1, p.y);
    Point right = Point(p.x+1, p.y);
    Point above_left  = Point(p.x-1, p.y+1);
    Point below_left  = Point(p.x-1, p.y-1);
    Point above_right = Point(p.x+1, p.y+1);
    Point below_right = Point(p.x+1, p.y-1);

    double phix = -phi.at<double>(left) + phi.at<double>(right);
    double phiy = -phi.at<double>(below) + phi.at<double>(above);
    double phixx = phi.at<double>(left) - 2*phi.at<double>(p) + phi.at<double>(right);
    double phiyy = phi.at<double>(below) - 2*phi.at<double>(p) + phi.at<double>(above);
    double phixy = -0.25*phi.at<double>(below_left) - 0.25*phi.at<double>(above_right)
                    +0.25*phi.at<double>(below_right) +0.25*phi.at<double>(above_left);
    double phix2 = pow(phix, 2);
    double phiy2 = pow(phiy, 2);

    //cout << phix << " " << phiy << " " << phixx << " " << phiyy << " " << phixy << " " << phix2 << " " << phiy2 << endl;
    return ((phix2*phiyy + phiy2*phixx - 2*phix*phiy*phixy) /
            pow(phix2 + phiy2 + eps, 1.5)) * pow((phix2 + phiy2), 0.5);
}

//
// Computes a single iteration of the Lankton active contour algorithm.
// Returns true if the contour converged, false otherwise.
//
bool activeContour(Mat src, Mat& phi) {

    double threshold_tolerance = 10e-16;  // stopping condition

    int radius = 10;             // radius of localized region to consider
    double alpha = 0.0;          // weight of local image stats energy
    double lambda = 1.0;         // weight of curvature energy
    double eps = 1.0;

    // get the curves narrow band
    vector<Point> border_points = getBorderPoints(phi);

    double total_energy = 0.0;
    double phi_diff = 0.0;
    // compute local energies along border of contour
    vector<Point> local_points;
    vector<double> image_energy;
    vector<double> curve_energy;

    double local_image_energy = 0.0;
    double ux, vx = 0.0;
    int I;
    for (auto p: border_points) {
        local_points = localSurroundingRegion(src, p, radius);
        ux = localRegionAverage(src, phi, local_points, INSIDE);
        vx = localRegionAverage(src, phi, local_points, OUTSIDE);

        // compute local stats and get image-based forces
        local_image_energy = 0.0;
        for (auto py: local_points) {
            // (I(y) - u_x)^2 - (I(y) - v_x)^2
            I = src.at<uchar>(py);
            local_image_energy += -(ux - vx) * (2*I - ux - vx);
        }
        image_energy.push_back(local_image_energy);

        curve_energy.push_back(getCurvature(phi, p, eps));
    }

    /*
    double max_energy = *max_element(begin(energies), end(energies));
    for (int i = 0; i < border_points.size(); ++i ) {
        Point p = border_points[i];
        total_energy = energies[i] / max_energy;
        cout << total_energy << endl;
        phi_diff += abs(phi.at<double>(p) - total_energy);
        phi.at<double>(p) = total_energy;
    }
    */

    vector<double> dphidt;

    double energy = 0.0;
    double max_image_energy = *max_element(image_energy.begin(), image_energy.end());
    double max_curve_energy = *max_element(curve_energy.begin(), curve_energy.end());
    max_image_energy = 1.0;
    //max_curve_energy = 1.0;
    for (int i = 0; i < border_points.size(); ++i ) {
        energy = alpha*(image_energy[i]/max_image_energy)
                + lambda*(curve_energy[i]/max_curve_energy);
        cout << curve_energy[i] << " " << energy << endl;
        dphidt.push_back(energy);
    }


    double max_dphidt = *max_element(dphidt.begin(), dphidt.end());
    // update phi and calculate total difference in energy
    double dt = 0.45 / (max_dphidt + eps);
    for (int i = 0; i < border_points.size(); ++i ) {
        Point p = border_points[i];
        total_energy = phi.at<double>(p) + (dphidt[i]);//dt*dphidt[i];
        //cout << "total energy " << total_energy << endl;
        phi_diff += abs(phi.at<double>(p) - total_energy);
        phi.at<double>(p) = total_energy;
    }

    // normalize total difference in old vs new phi (describes the contour)
    phi_diff = phi_diff / (phi.rows * phi.cols);
    cout << "difference: " << phi_diff << endl;

    // optionally reinitialize phi

    return (phi_diff < threshold_tolerance);
}

void getContourMaskFromPhi(Mat phi, Mat& contourMask) {
    contourMask = Mat::zeros(phi.size(), CV_8UC1);
    vector<Point> contour = foreground<double>(phi);
    for (auto point: contour) {
        contourMask.at<uchar>(point.y, point.x) = 255;
    }
}

