//opencv libraries
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

using namespace cv;
using namespace std;

#include "analysis.hpp"

vector<Point> read_landmarks_from_file(string filename) {
    vector<Point> points;

    ifstream f;
    f.open(filename);
    int x, y;
    while (f >> x >> y) {
        // create new point and add as result
        Point datapoint = Point(x, y);
        points.push_back(datapoint);
    }

    return points;
}

void rotMat(Mat& m, 

Point operator*(Mat M, const Point& p) {
    cv::Mat<double> src(3, 1);

    src(0,0) = p.x;
    src(1,0) = p.y;
    src(2,0) = 1.0;

    Mat dst = M*src;
    return Point(dst(0, 0), dst(1, 0));
}

void draw_contour(Mat& dst, vector<Point> points, Mat rot, Mat trans) {
    Point p1;
    Point p2;
    Mat result;
    for (int i = 1; i < points.size(); i++ ) {
        p1 = points[i-1];
        p2 = points[i];
		//p1 = rot.empty() ? p1, p1*rot;
		//p2 = rot.empty() ? p1, p1*rot;
		//p1 = trans.empty() ? p1, p1*trans;
		//p2 = trans.empty() ? p1, p1*trans;
        line(dst, p1, p2, Scalar(255), 2);
    }
    line(dst, points[points.size() - 1], points[0], Scalar(255), 2);
}

int main(int argc, char** argv)
{
    Mat img1;
    Mat img2;
    vector<Point> landmarks1;
    vector<Point> landmarks2;
    if (argc < 3) {
        cout << "Using default data files." << endl;
        landmarks1 = read_landmarks_from_file("../data/orientation1.txt");
        landmarks2 = read_landmarks_from_file("../data/orientation2.txt");
    } else {
        landmarks1 = read_landmarks_from_file(argv[1]);
        landmarks2 = read_landmarks_from_file(argv[2]);
    }
    assert(landmarks1.size() == landmarks2.size());

    // create images of points
    img1 = Mat::ones(400, 400, CV_8UC1);
    img2 = Mat::ones(400, 400, CV_8UC1);
    draw_contour(img1, landmarks1);
    draw_contour(img2, landmarks2);

    // merge into single vector, used to compute total center
    vector<Point> all_points;
    all_points.reserve(landmarks1.size() + landmarks2.size());
    all_points.insert(all_points.end(), landmarks1.begin(), landmarks1.end());
    all_points.insert(all_points.end(), landmarks2.begin(), landmarks2.end());

    // get center points
    Point lcenter = centroid(landmarks1);
    Point rcenter = centroid(landmarks2);
    Point center = centroid(all_points);

    Mat trans_mat = (Mat_<double>(2,3) << 1, 0, center.x, 0, 1, center.y);
    warpAffine(img1, img1, trans_mat, img1.size());
    warpAffine(img2, img2, trans_mat, img1.size());

    /*
    // move all points to shared center point for rotation calculation
    // this allows for rotation to be the only unknown in our equation

    // R = rotation for alignment
    // ro = translation for alignment
    // ro = rr_bar - R*rl_bar

    // get rotation R for alignment
    int num = 0;
    int denom = 0;
    int xl, yl, xr, yr;
    for (int i = 0; i < landmarks1.size(); ++i) {
        xl = center.x + (landmarks1[i].x - lcenter.x);
        yl = center.y + (landmarks1[i].y - lcenter.y);
        xr = center.x + (landmarks2[i].x - rcenter.x);
        yr = center.y + (landmarks2[i].y - rcenter.y);
        num += (yr*xl - xr*yl);
        denom += (xr*xl - yr*yl);
    }
    double theta = atan(num/denom);
    // R = |cos(theta)   -sin(theta)|
    //     |sin(tehta)    cos(theta)|
    double cos_r = cos(theta);
    double sin_r = sin(theta);

    // apply rotation, R*rl_bar
    int rot_lx = lcenter.x*cos_r - lcenter.y*sin_r;
    int rot_ly = lcenter.x*sin_r + lcenter.y*cos_r;

    // solve for translation, ro
    int ro_x = rcenter.x - rot_lx;
    int ro_y = rcenter.y - rot_ly;
    Point ro = Point(ro_x, ro_y);

    Mat translated = Mat::zeros(img1.size(), img1.type());
    Mat trans_mat = (Mat_<double>(2,3) << 1, 0, ro_x, 0, 1, ro_y);
    warpAffine(img1, translated, trans_mat, translated.size());

    Mat rotated = Mat::zeros(img1.size(), img1.type());
    Mat rot_mat = getRotationMatrix2D(rcenter, theta*(180/M_PI), 1.0);
    warpAffine(translated, rotated, rot_mat, rotated.size());
    */

    // Apply rotation and translation to first image to see
    // that it worked
    putText(img1, "image1", Point(10, img1.rows-10), 0, 0.5, Scalar(0,0,0));
    putText(img2, "image2", Point(10, img1.rows-10), 0, 0.5, Scalar(0,0,0));
    //putText(rotated, "image1 aligned to image2", Point(10, img1.rows-10), 0, 0.5, Scalar(0,0,0));
    Mat output = Mat(img1.rows, img1.cols * 3, CV_8UC1);
    img1.copyTo(output(Rect(0, 0, img1.cols, img1.rows)));
    img2.copyTo(output(Rect(img1.cols, 0, img1.cols, img1.rows)));
    //translated.copyTo(output(Rect(img1.cols*2, 0, img1.cols, img1.rows)));

    imshow("Output", output);
    waitKey(0);

    return 0;
}

