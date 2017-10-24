#ifndef __edgeac_H_INCLUDED__
#define __edgeac_H_INCLUDED__

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

double gradientEnergy(Mat& src, Point p);
double continuityEnergy(Point prev, Point p, double avg_dist);
double curvatureEnergy(Point l, Point m, Point n);
double curvatureFeedback(Point l, Point  m, Point n);
bool activeContour(Mat& src, vector<Point>& contour,
                    vector<double>& alpha,   // continuity weight
                    vector<double>& beta,    // curvature weight
                    vector<double>& gamma);  // gradient weight

#endif
