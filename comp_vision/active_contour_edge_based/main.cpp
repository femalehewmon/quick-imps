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
#include "edge_based_active_contour.hpp"

using namespace cv;
using namespace std;

// 1992 Williams and Shah active contour implementation

void drawActiveContour(Mat& dst, vector<Point> points, Mat edges) {
    if (points.size() <= 0) return;

    for (int i = 0; i < points.size(); ++i) {
        line(dst, points[i], points[(i+1) % points.size()], Scalar(0, 0, 255), 1);
        circle(dst, points[i], 2, Scalar(255,0,0), -1);
    }
    //circle(dst, points[0], 10, Scalar(100), -1);
    //circle(dst, points[points.size() - 1], 10, Scalar(100));
}

void drawMessage(Mat& dst, String msg) {
    putText(dst, msg, Point(20, 20), FONT_HERSHEY_SIMPLEX, 0.5, 100);
}

int main(int argc, char** argv)
{
    String img_file, init_file, edge_file;
    if (argc < 3) {
        img_file = "../data/image.png";
        init_file = "../data/init_doubled.txt";
    } else {
        img_file = argv[1];
        init_file = argv[2];
        edge_file = argv[3];
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
    //Canny(img, edges, 5, 150);
    edges = imread(edge_file, IMREAD_GRAYSCALE);

    // initialize weights for energy functions
    vector<double> alpha(contour.size(), 1.0); // continuity weights
    vector<double> beta(contour.size(), 1.0);  // curvature weights
    vector<double> gamma(contour.size(), 1.0); // gradient weights

    Mat composite;// = img.clone();
    cvtColor(img.clone(), composite, CV_GRAY2BGR);

    Mat buffer;
    cvtColor(edges.clone(), buffer, CV_GRAY2BGR);
    threshold(buffer, buffer, 1, 255, THRESH_BINARY);
    buffer.copyTo(composite, buffer);

    imwrite("orig.png", composite);
    Mat vis = composite.clone();
    drawActiveContour(vis, contour, edges);
    drawMessage(vis, "Waiting to start... press any key");

    imwrite("start.png", vis);
    imshow("Output", vis);
    waitKey(0);

    // run active contour algorithm and visualized progress
    int it_count = 0;
    while (it_count < 200 && activeContour(edges, contour, alpha, beta, gamma)) {
        // increase iteration count
        ++it_count;
        cout << "ITERATION: " << it_count << endl;

        // visualize active contour iteration
        vis = composite.clone();
        drawActiveContour(vis, contour, edges);
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
    waitKey(0); // display final frame
}

