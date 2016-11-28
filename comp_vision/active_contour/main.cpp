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

using namespace cv;
using namespace std;

bool operator==(const Point& m, const Point& n) {
    return (m.x == n.x && m.y == n.y);
}
bool operator!=(const Point& m, const Point& n) {
    return !(m == n);
}

double gradient(Mat& src, Point p) {
    return 0.0;
}

double continuity(Mat& src, Point p) {
    return 0.0;
}

double curvature(Mat& src, Point p) {
    return 0.0;
}

bool activeContour(Mat& src, vector<Point>& contour) {
    double alpha = 1.0;
    double beta  = 1.0;
    double gamma = 1.0;

    int threshold = contour.size(); // stop after 1 iteration

    int points_moved = 0;
    int n = contour.size();

    int i, j;
    double Ej, Emin, Econt, Ecurv, Eimage = 0.0;
    int jmin;
    vector<Point> neighbors;
    for ( i = 0; i < n; ++i ) {
        // for each point, determine if a location in the neighborhood
        // would minimize the active contour energies
        jmin = 0;
        Emin = numeric_limits<double>::max();
        neighbors = neighborhood(contour[i], N4, true);
        for ( j = 0; j < neighbors.size(); ++j ) {
            Econt = continuity(src, neighbors[j]);
            Ecurv = curvature(src, neighbors[j]);
            Eimage = gradient(src, neighbors[j]);
            Ej = alpha*Econt + beta*Ecurv + gamma*Eimage;
            if (Ej < Emin) {
                Emin = Ej;
                jmin = j;
            }
        }
        // if a better location is found for the current point, move it
        if ( neighbors[jmin] != contour[i] ) {
            contour[i] = neighbors[jmin];
            ++points_moved;
        }
    }
    // determine where to allow corners in the next iteration

    cout << points_moved << endl;
    return (points_moved >= threshold);
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

    Mat img = imread(img_file, IMREAD_GRAYSCALE);
    if ( !img.data ) {
        cout << "Could not open file " << img_file << endl;
        return -1;
    }
    vector<Point> contour = loadPointsFromFile(init_file);
    if ( contour.size() <= 0 ) {
        cout << "Could not open file " << init_file << endl;
        return -1;
    }

    drawVector(img, contour);

    bool iterate = true;
    while (iterate) {
        iterate = activeContour(img, contour);
    }

    imshow("Output", img);
    waitKey(0);
}

