#ifndef __utils_H_INCLUDED__
#define __utils_H_INCLUDED__

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

vector<Point> loadPointsFromFile(string filename);
void drawVector(Mat& dst, vector<Point> points);

enum NEIGHBOR {N4, N8};
vector<Point> neighborhood(Point p, NEIGHBOR type, bool include_self=false);

#endif
