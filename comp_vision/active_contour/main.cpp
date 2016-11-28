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

using namespace cv;
using namespace std;

// William and Shah active contour implementation

bool operator==(const Point& m, const Point& n) {
    return (m.x == n.x && m.y == n.y);
}
bool operator!=(const Point& m, const Point& n) {
    return !(m == n);
}

double gradientEnergy(Mat& src, Point p) {
    // first order derivative in x direction
    double di_dx = 0;
    di_dx += src.at<uchar>(p.y, p.x+1) - src.at<uchar>(p.y, p.x);
    di_dx += src.at<uchar>(p.y+1, p.x+1) - src.at<uchar>(p.y+1, p.x);
    di_dx *= 0.5;
    // first order derivative in y direction
    double di_dy = 0;
    di_dy += src.at<uchar>(p.y+1, p.x) - src.at<uchar>(p.y, p.x);
    di_dy += src.at<uchar>(p.y+1, p.x+1) - src.at<uchar>(p.y, p.x+1);
    di_dy *= 0.5;
    cout << "change in x, y " << di_dx << " " << di_dy << endl;
    // squared magnitude of gradient
    return pow(sqrt(pow(di_dx, 2) + pow(di_dy, 2)), 2);
}

double continuityEnergy(Point prev, Point p, double avg_dist) {
    // unnormalized value
    double curr_dist = distanceBetweenPoints(prev, p);
    return abs(avg_dist - curr_dist);
}


double curvatureEnergy(Point l, Point m, Point n) {
    // unnormalized value
    curvature(l, m, n);
}

bool activeContour(Mat& src, vector<Point>& contour) {
    double alpha = 1.0;
    double beta  = 1.0;
    double gamma = 1.0;

    int threshold = 1; // stop after 1 iteration

    int points_moved = 0;
    int n = contour.size();
    vector<double> distances = distanceBetweenPoints(contour);
    double avg_dist = accumulate(distances.begin(), distances.end(), 0)/n;
    double max_dist = 0.0;

    int i, j;
    double Ej, Emin, Econt, Ecurv, Eimage = 0.0;
    int jmin;
    double gmin, gmax;
    int mag;
    vector<Point> neighbors;
    for ( i = 1; i < n-1; ++i ) {
        // for each point, determine if a location in the neighborhood
        // would minimize the active contour energies
        jmin = 0;
        Emin = numeric_limits<double>::max();
        neighbors = neighborhood(contour[i], N8, true);

        // calculate continuity and curvature energy of each neighbor
        vector<double> continuity;
        vector<double> curvature;
        vector<double> gradient;
        for ( j = 0; j < neighbors.size(); ++j ) {
            // continuity
            Econt = continuityEnergy(neighbors[j], contour[i-1], avg_dist);
            continuity.push_back(Econt);
            // curvature
            Ecurv = curvatureEnergy(contour[i-1], neighbors[j], contour[i+1]);
            curvature.push_back(Ecurv);
            // gradient
            Eimage = gradientEnergy(src, neighbors[j]);
            gradient.push_back(Eimage);
        }
        // normalize values to be used in full energy function
        normalize(continuity);
        normalize(curvature);
        minMaxLoc(gradient, &gmin, &gmax);

        // calculate final energy function and move points if needed
        for ( j = 0; j < neighbors.size(); ++j ) {
            Econt = continuity[j];
            Ecurv = curvature[j];
            // normalize gradient energy in place
            mag = src.at<uchar>(neighbors[j].y, neighbors[j].y);
            Eimage = (gmax - gmin) > 0 ? (gmin - mag)/(gmax - gmin) : 0;
            cout << "image mag: " << mag << endl;
            cout << "max and min grad: " << gmax << " " << gmin << endl;

            cout << "continuity: " << Econt << endl;
            cout << "curvature: " << Ecurv << endl;
            cout << "gradient: " << Eimage << endl;
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
    // add start point to end of contour and vise versa to complete loop
    contour.insert(contour.begin(), contour[contour.size() - 1]);
    contour.push_back(contour[1]);

    int it_max = 200;

    Mat vis;
    int it_count = 0;
    while (activeContour(img, contour)) {
        // visualize active contour iteration
        vis = img.clone();
        drawVector(vis, contour);
        imshow("Output", vis);
        waitKey(30);

        // increase iteration count
        ++it_count;
        if (it_count >= it_max) break;
    }

}

