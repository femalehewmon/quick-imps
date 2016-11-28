#ifndef __vectorize_H_INCLUDED__
#define __vectorize_H_INCLUDED__

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

using namespace cv;
using namespace std;

vector<Point> binaryImageToVectorizedContour(Mat& src, int obj_id=255);
vector<Point> floodFillAndVectorize(Mat& src, Mat& buffer,
                                    vector<Point> searchSpace, int found=255);
void neighborhood(vector<Point> points, vector<Point>& neighbors);
bool inBounds(Point p, int w, int h);
Point borderPoint(Mat& src, int obj_id=255);

#endif
