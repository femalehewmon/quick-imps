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
    if ( points.size() <= 0 ) return;
    points.push_back(points[0]);
    for (int i = 0; i < points.size() - 1; i++ ) {
        line(dst, points[i], points[i+1], Scalar(255), 2);
    }
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

