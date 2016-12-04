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

Point operator*(Mat M, const Point& p) {
    Mat src = (Mat_<uchar>(3, 1) << p.x, p.y, 1.0);
    Mat dst = M*src;
    return Point(dst.at<uchar>(0, 0), dst.at<uchar>(1, 0));
}

void draw_contour(Mat& dst, vector<Point> points, Mat trans=Mat(), Mat rot=Mat());
void draw_contour(Mat& dst, vector<Point> points, Mat trans, Mat rot) {
    Point p1;
    Point p2;
    Mat result;
    points.push_back(points[0]);
    for (int i = 1; i < points.size(); i++ ) {
        p1 = points[i-1];
        p2 = points[i];
        /*
		p1 = rot.empty() ? p1 : rot*p1;
		p2 = rot.empty() ? p2 : rot*p2;
		p1 = trans.empty() ? p1 : trans*p1;
		p2 = trans.empty() ? p1 : trans*p2;
        */
        line(dst, p1, p2, Scalar(255), 2);
    }
}
void draw_contour2(Mat& dst, vector<Point> points) {
    Point p1, p2;
    points.push_back(points[0]);
    for (int i = 1; i < points.size(); i++ ) {
        p1 = points[i-1];
        p2 = points[i];
        line(dst, p1, p2, Scalar(100), 2);
    }
}

void mergePoints(vector<Point>& all_points, vector<Point> p1, vector<Point> p2) {
    all_points.reserve(p1.size() + p2.size());
    all_points.insert(all_points.end(), p1.begin(), p1.end());
    all_points.insert(all_points.end(), p2.begin(), p2.end());
}

int main(int argc, char** argv)
{
    Mat img1;
    vector<Point> landmarks;
    if (argc < 3) {
        cout << "Using default data files." << endl;
        landmarks = read_landmarks_from_file("../data/orientation1.txt");
    } else {
        landmarks = read_landmarks_from_file(argv[2]);
    }
    /*
    assert(landmarksLeft.size() == landmarksRight.size());

    // shift second set of landmarks
    for (int i = 0; i < landmarksRight.size(); ++i ) {
        landmarksRight[i] = Point(landmarksRight[i].x + 400, landmarksRight[i].y);
    }
    */
    int x, y;

    vector<Point> a = landmarks;
    vector<Point> b = landmarks;
    Point center = centroid(a);
    for (int i = 0; i < b.size(); ++i) {
        x = 0 + (b[i].x - center.x);
        y = 0 + (b[i].y - center.y);
        b[i] = Point(x, y);
    }
    // rotate b
    Mat rotMat = getRotationMatrix2D(center, 45, 1.0);
    for (int i = 0; i < b.size(); ++i) {
        x = rotMat.at<double>(0,0)*b[i].x + rotMat.at<double>(0,1)*b[i].y;
        y = rotMat.at<double>(1,0)*b[i].x + rotMat.at<double>(1,1)*b[i].y;
        b[i] = Point(x, y);
    }
    for (int i = 0; i < b.size(); ++i) {
        x = center.x + (b[i].x);
        y = center.y + (b[i].y);
        b[i] = Point(x, y);
    }
    Mat output = Mat(400, 400, CV_8UC1);
    draw_contour(output, a);
    draw_contour(output, b);

    // apply rotation, R*rl_bar
    int rot_lx = lcenter.x*cos_r - lcenter.y*sin_r;
    int rot_ly = lcenter.x*sin_r + lcenter.y*cos_r;

    // solve for translation, ro
    int ro_x = rcenter.x - rot_lx;
    int ro_y = rcenter.y - rot_ly;
    Point ro = Point(ro_x, ro_y);

    // get rotation R for alignment
    int num = 0;
    int denom = 0;
    int xl, yl, xr, yr;
    for (int i = 0; i < b.size(); ++i) {
        num += (b[i].y*a[i].x - b[i].x*a[i].y);
        denom += (b[i].x*a[i].x - b[i].y*a[i].y);
    }
    double theta = atan((double)num/(double)denom)*(180/M_PI);
    cout << theta << endl;

    rotMat = getRotationMatrix2D(center, theta, 1.0);
    for (int i = 0; i < b.size(); ++i) {
        x = 0 + (b[i].x - center.x);
        y = 0 + (b[i].y - center.y);
        b[i] = Point(x, y);
    }
    for (int i = 0; i < b.size(); ++i) {
        x = rotMat.at<double>(0,0)*b[i].x + rotMat.at<double>(0,1)*b[i].y;
        y = rotMat.at<double>(1,0)*b[i].x + rotMat.at<double>(1,1)*b[i].y;
        b[i] = Point(x, y);
    }
    for (int i = 0; i < b.size(); ++i) {
        x = center.x + (b[i].x);
        y = center.y + (b[i].y);
        b[i] = Point(x, y);
    }
    draw_contour2(output, b);

    imshow("Output", output);
    waitKey(0);
    return 0;


    /*

    vector<Point> all_points;
    mergePoints(all_points, landmarksLeft, landmarksRight);

    // create images of points
    img1 = Mat::ones(400, 400, CV_8UC1);
    img2 = Mat::ones(400, 400, CV_8UC1);
    Mat output = Mat(img1.rows, img1.cols * 2, CV_8UC1);

    draw_contour(output, landmarksLeft);
    draw_contour(output, landmarksRight);

    // get center points
    Point lcenter = centroid(landmarksLeft);
    Point rcenter = centroid(landmarksRight);
    Point center = centroid(all_points);            // shared center

    drawMarker(output, lcenter, Scalar(100), MARKER_CROSS);
    drawMarker(output, rcenter, Scalar(100), MARKER_CROSS);
    drawMarker(output, center, Scalar(100), MARKER_CROSS);

    // move all points to shared center point for rotation calculation
    // this allows for rotation to be the only unknown in our equation

    // R = rotation for alignment
    // ro = translation for alignment
    // ro = rr_bar - R*rl_bar

    // get rotation R for alignment
    int num = 0;
    int denom = 0;
    int xl, yl, xr, yr;
    for (int i = 0; i < landmarksLeft.size(); ++i) {
        xl = center.x + (landmarksLeft[i].x - lcenter.x);
        yl = center.y + (landmarksLeft[i].y - lcenter.y);
        landmarksLeft[i] = Point(xl, yl);
        xr = center.x + (landmarksRight[i].x - rcenter.x);
        yr = center.y + (landmarksRight[i].y - rcenter.y);
        landmarksRight[i] = Point(xr, yr);
        num += (yr*xl - xr*yl);
        denom += (xr*xl - yr*yl);
    }
    for (int i = 0; i < landmarksLeft.size(); ++i) {
        xl = 0 + (landmarksLeft[i].x - center.x);
        yl = 0 + (landmarksLeft[i].y - center.y);
        landmarksLeft[i] = Point(xl, yl);
    }
    double theta = atan((double)num/(double)denom);//*(180/M_PI);
    cout << theta << endl;
    // R = |cos(theta)   -sin(theta)|
    //     |sin(tehta)    cos(theta)|
    double cos_r = cos(theta);
    double sin_r = sin(theta);

    Point p;
    double x, y;
    //Mat rotMat = getRotationMatrix2D(center, theta, 1.0);
    Mat rotMat = (Mat_<double>(2, 2) << cos_r, -sin_r, sin_r, cos_r);
    for (int i = 0; i < landmarksLeft.size(); ++i) {
        p = landmarksLeft[i];
        x = rotMat.at<double>(0,0)*p.x + rotMat.at<double>(0,1)*p.y;
        y = rotMat.at<double>(1,0)*p.x + rotMat.at<double>(1,1)*p.y;
        landmarksLeft[i] = Point(x, y);
    }
    draw_contour2(output, landmarksLeft);
    //draw_contour2(output, landmarksRight);
    for (int i = 0; i < landmarksLeft.size(); ++i) {
        xl = center.x + (landmarksLeft[i].x);
        yl = center.y + (landmarksLeft[i].y);
        landmarksLeft[i] = Point(xl, yl);
    }
    draw_contour2(output, landmarksLeft);

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

    imshow("Output", output);
    waitKey(0);
    return 0;
    */
}

