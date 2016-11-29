//opencv libraries
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/core/core.hpp"
//C++ standard libraries
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <glob.h>
#include <math.h>
#include <limits>

#include "utils.hpp"
#include "characterization.hpp"

using namespace cv;
using namespace std;

// 1992 Williams and Shah active contour implementation
//
// Williams, Donna J., and Mubarak Shah.
//     "A fast algorithm for active contours and curvature estimation."
//     CVGIP: Image understanding 55.1 (1992): 14-26.

double gradientEnergy(Mat& src, Point p) {
    /*
    // unnormalized value
    // first order derivative in x direction
    double di_dx = 0;
    di_dx += src.at<uchar>(p.y, p.x+1) - src.at<uchar>(p.y, p.x);
    di_dx += src.at<uchar>(p.y+1, p.x+1) - src.at<uchar>(p.y+1, p.x);
    di_dx *= 0.5;
    // first order derivative in y direction
    double di_dy = 0;
    di_dy += src.at<uchar>(p.y+1, p.x) - src.at<uchar>(p.y, p.x);
    di_dy += src.at<uchar>(p.y+1, p.x+1) - src.at<uchar>(p.y, p.x+1);
    di_dy *= 0.5;
    // squared magnitude of gradient
    return pow(sqrt(pow(di_dx, 2) + pow(di_dy, 2)), 2);
    */
    // return magnitude only since edges are passed in as input
    return src.at<uchar>(p.y, p.x);
}

double continuityEnergy(Point prev, Point p, double avg_dist) {
    // unnormalized value
    double curr_dist = distanceBetweenPoints(prev, p);
    return abs(avg_dist - curr_dist);
}

double curvatureEnergy(Point l, Point m, Point n) {
    // unnormalized value
    // |vss|^2
    return pow(l.x - 2*m.x + n.x, 2) + pow(l.y - 2*m.y + n.y, 2);
}

double curvatureFeedback(Point l, Point  m, Point n) {
    double dxi = (m.x - l.x);
    double dxi1 = (n.x - m.x);
    double dyi = (m.y - l.y);
    double dyi1 = (n.y - m.y);
    double dsi = sqrt(pow(dxi,2) + pow(dyi,2));
    double dsi1 = sqrt(pow(dxi1,2) + pow(dyi1,2));

    return pow((dxi/dsi) - (dxi1/dsi1), 2) + pow((dyi/dsi) - (dyi1/dsi1), 2);
}

bool activeContour(Mat& src, vector<Point>& contour,
                    vector<double>& alpha,    // continuity weight
                    vector<double>& beta,     // curvature weight
                    vector<double>& gamma ){  // gradient weight

    int n = contour.size(); // total number of points
    int points_moved = 0;   // points moved this iteration

    int threshold_moved = 1; // minimum number of points that must
                             // update to continue iterating

    // continuity energy calculation helpers
    vector<double> distances = distanceBetweenPoints(contour);
    double avg_dist = accumulate(distances.begin(), distances.end(), 0) / n;
    double max_dist = 0.0;

    int i, j, jmin = 0;
    double Ej, Emin, Econt, Ecurv, Eimage = 0.0;
    double gmin, gmax = 0.0;
    int mag = 0;
    int prev, curr, next = 0;
    vector<Point> neighbors;
    // for each point, determine if a location in its neighborhood
    // would minimize the active contour energies
    for (i = 0; i <= n; ++i) { // process first point twice
        prev = (i-1) % n;
        curr = i % n;
        next = (i+1) % n;

        jmin = 0;
        Emin = numeric_limits<double>::max();
        neighbors = neighborhood(contour[curr], N8, true);

        // calculate un-normalized energies for each potential neighbor point
        vector<double> continuity;
        vector<double> curvature;
        vector<double> gradient;
        for (j = 0; j < neighbors.size(); ++j) {
            // continuity
            Econt = continuityEnergy(neighbors[j], contour[prev], avg_dist);
            continuity.push_back(Econt);
            // curvature
            Ecurv = curvatureEnergy(contour[prev], neighbors[j], contour[next]);
            curvature.push_back(Ecurv);
            // gradient
            Eimage = gradientEnergy(src, neighbors[j]);
            gradient.push_back(Eimage);
        }
        // normalize values to be used in full energy function
        normalize(continuity);
        normalize(curvature);
        minMaxLoc(gradient, &gmin, &gmax); // helpers to normalize gradient

        // calculate final energy function and move points if needed
        for (j = 0; j < neighbors.size(); ++j) {
            Econt = continuity[j];
            Ecurv = curvature[j];
            // normalize gradient energy in place
            mag = src.at<uchar>(neighbors[j].y, neighbors[j].x);
            Eimage = (gmax - gmin) > 0 ? (gmin - mag)/(gmax - gmin) : 0;

            //cout << i << " neighbor " << j << endl;
            //cout << "continuity: " << Econt << endl;
            //cout << "curvature: " << Ecurv << endl;
            //cout << "gradient: " << Eimage << endl;

            // total energy function
            Ej = alpha[curr]*Econt + beta[curr]*Ecurv + gamma[curr]*Eimage;

            // update min point if energy is minimized
            if (Ej < Emin) {
                Emin = Ej;
                jmin = j;
            }
        }
        // if a better location was found for the current point, move it
        if (neighbors[jmin] != contour[curr]) {
            contour[curr] = neighbors[jmin];
            ++points_moved;
        }
    }

    // determine where to allow corners in the next iteration
    vector<double> cfeedback;
    for (i = 0; i < n; ++i) {
        prev = (i-1) % n;
        next = (i+1) % n;
        cfeedback.push_back(
                curvature(contour[prev], contour[i], contour[next]));
    }
    double threshold_cfeedback = 0;
    double threshold_magnitude = 100;
    for (i = 0; i < n; ++i) {
        prev = (i-1) % n;
        next = (i+1) % n;
        mag = src.at<uchar>(contour[i].y, contour[i].x);
        if (cfeedback[i] > cfeedback[prev] && cfeedback[i] > cfeedback[next]
                && cfeedback[i] > threshold_cfeedback
                && mag > threshold_magnitude) {
            beta[i] = 0;
        }
    }

    cout << "points moved: " << points_moved << " out of " << n << endl;
    return (points_moved >= threshold_moved);
}

void drawActiveContour(Mat& dst, vector<Point> points) {
    if (points.size() <= 0) return;
    circle(dst, points[0], 3, Scalar(0), -1);
    for (int i = 1; i < points.size(); ++i) {
        circle(dst, points[i], 3, Scalar(100));
        //line(dst, points[i], points[i+1], Scalar(100), 2);
    }
    circle(dst, points[points.size() - 1], 10, Scalar(100));
}

void drawMessage(Mat& dst, String msg) {
    putText(dst, msg, Point(20, 20), FONT_HERSHEY_SIMPLEX, 0.5, 100);
}

int main(int argc, char** argv)
{
    String img_file, init_file;
    if (argc < 3) {
        img_file = "../data/image.png";
        init_file = "../data/init.txt";
    } else {
        img_file = argv[1];
        init_file = argv[2];
    }

    Mat img = imread(img_file, IMREAD_GRAYSCALE);
    if ( !img.data ) {
        cout << "Could not open file " << img_file << endl;
        return -1;
    }
    vector<Point> contour = loadPointsFromFile(init_file);
    if ( contour.size() <= 0 ) {
        cout << "Could not read points from file " << init_file << endl;
        return -1;
    }

    // run edge detection as input to active contour
    Mat edges;
    Canny(img, edges, 0, 100);

    // initialize weights for energy functions
    vector<double> alpha(contour.size(), 1.0); // continuity weights
    vector<double> beta(contour.size(), 1.0);  // curvature weights
    vector<double> gamma(contour.size(), 1.0); // gradient weights

    Mat vis = img.clone();
    drawActiveContour(vis, contour);
    drawMessage(vis, "Waiting to start... press any key");
    imshow("Output", vis);
    waitKey(0);

    // run active contour algorithm and visualized progress
    int it_count = 0;
    while (it_count < 150 && activeContour(edges, contour, alpha, beta, gamma)) {
        // increase iteration count
        ++it_count;
        cout << "ITERATION: " << it_count << endl;

        // visualize active contour iteration
        vis = img.clone();
        drawActiveContour(vis, contour);
        drawMessage(vis, "iteration: " + to_string(it_count));
        imshow("Output", vis);
        waitKey(30);
    }
    waitKey(0); // display final frame

}

