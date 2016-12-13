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
    Mat contourMask;
    getContourMaskFromPhi(phi, contourMask);

	Mat edges;
    vector<vector<Point> > contours;
    vector<Vec4i> hierarchy;
    Canny(contourMask, edges, 1, 100, 3);
    findContours(edges, contours, hierarchy,
                        CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE);
    for( int i = 0; i< contours.size(); i++ ) {
        drawContours(dst, contours, i, Scalar(0, 0, 255), 2, 8, hierarchy, 0);
    }
    /*
    vector<Point> points = foreground<uchar>(contourMask);
    for (auto point: points) {
        circle(dst, point, 5, Scalar(255, 0, 0), -1);
    }
    */
}

void drawMessage(Mat& dst, String msg) {
    putText(dst, msg, Point(20, 20), FONT_HERSHEY_SIMPLEX, 0.5, 100);
}

int main(int argc, char** argv)
{
    String img_file;
    if (argc < 2) {
        img_file = "../data/image.png";
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

    Mat composite;
    cvtColor(img.clone(), composite, CV_GRAY2BGR);

    // show initial visualization, prompt user to begin segmentation
    Mat vis = composite.clone();
    drawActiveContour(vis, phi);
    drawMessage(vis, "Waiting to start... press any key");
    imwrite("start.png", vis);
    imshow("Output", vis);
    waitKey(0);

    // run active contour algorithm and visualize progress
    int it_max = 200; // cut-off point if algorithm does not converge
    int it_count = 0;
    while (it_count < it_max && !activeContour(img, phi)) {
        // increase iteration count
        ++it_count;
        cout << "ITERATION: " << it_count << endl;

        // visualize active contour iteration
        vis = composite.clone();
        drawActiveContour(vis, phi);
        drawMessage(vis, "iteration: " + to_string(it_count));
        if (it_count == 25) {
            imwrite("mid1.png", vis);
        } else if (it_count == 50) {
            imwrite("mid2.png", vis);
        } else if (it_count == 100) {
            imwrite("mid3.png", vis);
        }
        imshow("Output", vis);
        waitKey(30);
    }
    imwrite("end.png", vis);
    waitKey(0); // display final frame until key pressed
}

