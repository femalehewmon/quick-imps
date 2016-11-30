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

    // load the input image
    Mat img = imread(img_file, IMREAD_GRAYSCALE);
    if ( !img.data ) {
        cout << "Could not open file " << img_file << endl;
        return -1;
    }

    Mat init = imread(init_file, IMREAD_GRAYSCALE);
    if ( !img.data ) {
        cout << "Could not open file " << img_file << endl;
        return -1;
    }

    // initialize phi to binary input initialization
    // phi > 0 : inside contour
    // phi < 0 : outside contour
    // phi = 0 : on contour
    Mat phi = Mat::zeros(init.size(), CV_32F);
    threshold(init, phi, 1, 1, THRESH_BINARY);

    // initialize weights for energy functions
    //Mat phi = Mat(img.size(), CV_32F);
    //initializePhi(phi);     // initialize to checkerboard pattern

    // show initial visualization, prompt user to begin segmentation
    Mat vis = img.clone();
    drawActiveContour(vis, contour);
    drawMessage(vis, "Waiting to start... press any key");
    imshow("Output", vis);
    waitKey(0);

    // run active contour algorithm and visualize progress
    int it_max = 500; // cut-off point if algorithm does not converge
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

