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

vector<Point> read_datapoints_from_file(string filename) {
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

int main(int argc, char** argv)
{
    vector<Point> landmarks1;
    vector<Point> landmarks2;
    if (argc < 3) {
        cout << "Using default data files." << endl;
        landmarks1 = read_datapoints_from_file("../data/orientation1.txt");
        landmarks2 = read_datapoints_from_file("../data/orientation2.txt");
    } else {
        landmarks1 = read_datapoints_from_file(argv[1]);
        landmarks2 = read_datapoints_from_file(argv[2]);
    }

    return 0;
}

