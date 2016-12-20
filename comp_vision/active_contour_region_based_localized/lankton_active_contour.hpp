#ifndef __lankton_H_INCLUDED__
#define __lankton_H_INCLUDED__

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

void initializePhiCircle(Mat& phi, int width, int height);
double areaInsideContour(Mat phi);

enum RegionType {INSIDE, OUTSIDE};
double regionAverage(Mat src, Mat phi, vector<Point> B, double eps, RegionType region);

double heaviside (double phi_n, double eps);
double diracDelta (double phi_n, double eps);
double curveEnergy (Mat phi, Point p, double h, double mu);
bool activeContour(Mat src, Mat& phi);
void getContourMaskFromPhi(Mat phi, Mat& contourMask);


#endif
