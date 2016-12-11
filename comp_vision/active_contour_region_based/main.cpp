// opencv libraries
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/core/core.hpp"
// C++ standard libraries
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
// custom libraries
#include "utils.hpp"
#include "characterization.hpp"
#include "chan_vese_active_contour.hpp"

using namespace cv;
using namespace std;

//
// Demo of 2001 Chan-Vese active contour implementation
//
// Uses a built-in checker-board grid as the initialization of phi by default.
// If custom initialization curve desired, set flag when running.
//

void drawActiveContour(Mat& dst, Mat phi) {
    int row, col;
    for ( row = 0; row < phi.rows; ++row ) {
        for ( col = 0; col < phi.cols; ++col ) {
            if (abs(phi.at<double>(row, col)) == 0) {
                circle(dst, Point(col, row), 5, Scalar(0), -1);
            } else if (phi.at<double>(row, col) > 0 ) {
                circle(dst, Point(col, row), 1, Scalar(100));
            }
        }
    }
    return;
    /*
    if (points.size() <= 0) return;
    circle(dst, points[0], 3, Scalar(0), -1);
    for (int i = 1; i < points.size(); ++i) {
        circle(dst, points[i], 3, Scalar(100));
        //line(dst, points[i], points[i+1], Scalar(100), 2);
    }
    circle(dst, points[points.size() - 1], 10, Scalar(100));
    */
}

void drawMessage(Mat& dst, String msg) {
    putText(dst, msg, Point(20, 20), FONT_HERSHEY_SIMPLEX, 0.5, 100);
}

int main(int argc, char** argv)
{
    String img_file;
    if (argc < 2) {
        img_file = "../data/test1.jpg";
        //img_file = "../data/Snake-contour-example.jpg";
    } else {
        img_file = argv[1];
    }

    // load the input image
    Mat img = imread(img_file, IMREAD_GRAYSCALE);
    if ( !img.data ) {
        cout << "Could not open file " << img_file << endl;
        return -1;
    }
    cout << img.size() << endl;

    // initialize phi to initialization function
    // phi > 0 : inside contour
    // phi < 0 : outside contour
    // phi = 0 : on contour
    Mat phi = Mat::zeros(img.size(), CV_64F);
    initializePhiCircle(phi, phi.cols, phi.rows);
    //initializePhiCheckerboard(phi);

    // initialize weights for energy functions
    vector<Point> contour;
    getContourFromPhi(phi, contour);

    // show initial visualization, prompt user to begin segmentation
    Mat vis = img.clone();
    drawActiveContour(vis, phi);
    drawMessage(vis, "Waiting to start... press any key");
    imshow("Output", vis);
    waitKey(0);

    // run active contour algorithm and visualize progress
    int it_max = 100; // cut-off point if algorithm does not converge
    int it_count = 0;
    while (it_count < it_max && !activeContour(img, phi)) {
        // increase iteration count
        ++it_count;
        cout << "ITERATION: " << it_count << endl;

        // visualize active contour iteration
        vis = img.clone();
        drawActiveContour(vis, phi);
        drawMessage(vis, "iteration: " + to_string(it_count));
        imshow("Output", vis);
        waitKey(30);
    }
    waitKey(0); // display final frame until key pressed
}

