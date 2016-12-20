#include "utils.hpp"

vector<Point> loadPointsFromFile(string filename) {
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

void drawVector(Mat& dst, vector<Point> points) {
    if (points.size() <= 0) return;
    for (int i = 1; i < points.size() - 1; ++i) {
        line(dst, points[i], points[i+1], Scalar(100), 2);
    }
}

Point farthestPoint(Point p, vector<Point> ps) {
    if (ps.size() <= 0) return Point(-1, -1);
    int i, imax = 0;
    double d, dmax = 0.0;
    for (i = 0; i < ps.size(); ++i) {
        d = distanceBetweenPoints(p, ps[i]);
        if (d > dmax) {
            dmax = d;
            imax = i;
        }
    }
    return ps[i];
}

bool inBounds (Mat image, Point p) {
    return p.x >= 0 && p.y >= 0 && p.x < image.cols && p.y < image.rows;
}

double distanceBetweenPoints(Point a, Point b) {
    return sqrt(pow(b.x - a.x, 2) + pow(b.y - a.y, 2));
}

vector<double> distanceBetweenPoints(vector<Point> points) {
    // average distance between points
    vector<double> dists;
    if (points.size() <= 0) return dists;
    points.push_back(points[0]);

    for (int i = 0; i < points.size() - 1; ++i) {
        dists.push_back(distanceBetweenPoints(points[i], points[i+1]));
    }
    return dists;
}

vector<double> distanceBetweenPoints(vector<Point> points, Point ref_point) {
    // average distance between points
    vector<double> dists;
    if (points.size() <= 0) return dists;

    for (int i = 0; i < points.size() - 1; ++i) {
        dists.push_back(distanceBetweenPoints(ref_point, points[i]));
    }
    return dists;
}

vector<Point> neighborhood(Point p, NEIGHBOR type, bool include_self) {
    // stored in clockwise position starting at immediate left position
    // NOTE: if self include requested, it is added as the first point
    vector<Point> neighbors;
    if (include_self) neighbors.push_back(p);
    switch(type) {
        case N4:
            neighbors.push_back(Point(p.x-1, p.y));
            neighbors.push_back(Point(p.x, p.y-1));
            neighbors.push_back(Point(p.x+1, p.y));
            neighbors.push_back(Point(p.x, p.y+1));
            break;
        case N8:
            neighbors.push_back(Point(p.x-1, p.y));
            neighbors.push_back(Point(p.x-1, p.y-1));
            neighbors.push_back(Point(p.x, p.y-1));
            neighbors.push_back(Point(p.x+1, p.y-1));
            neighbors.push_back(Point(p.x+1, p.y));
            neighbors.push_back(Point(p.x+1, p.y+1));
            neighbors.push_back(Point(p.x, p.y+1));
            neighbors.push_back(Point(p.x-1, p.y+1));
            break;
    }
    return neighbors;
}

