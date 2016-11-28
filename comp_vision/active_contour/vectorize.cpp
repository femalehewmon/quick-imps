#include "characterization.hpp"
#include "vectorize.hpp"

vector<Point> binaryImageToVectorizedContour(Mat& src, int obj_id) {
    Mat buffer = Mat::zeros(src.size(), src.type());

    // dividing line is the initial search space on which to begin flood filling
    vector<Point> searchSpace;
    vector<Point> dividingLine;
    // create a dividing line through center of skeleton
    Point p;
    Point s = borderPoint(src, obj_id);

    int found = obj_id;
    vector<Point> leftSearch;
    vector<Point> rightSearch;

    vector<Point> obj = foreground(src, obj_id);
    double theta = orientation(obj);
    if (theta < 2.35 && theta > 0.785 ) {
        // assumes a vertical orientation
        for ( int col = s.x; col < src.cols - s.x; ++col ) {
            if ( src.at<uchar>(s.y, col) == obj_id ) {
                p = Point(col, s.y);
                dividingLine.push_back(p);
                searchSpace.push_back(p);
            } else {
                break;
            }
        }

        neighborhood(dividingLine, searchSpace);
        for (auto v = searchSpace.begin(); v != searchSpace.end(); ++v ) {
            p = *v;
            buffer.at<uchar>(p.y, p.x) = found;
            if ( p.y < s.y ) {
                leftSearch.push_back(p);
            } else if ( p.y > s.y ) {
                rightSearch.push_back(p);
            }
        }

    } else {
        // assumes a horizontal orientation
        for ( int row = s.y+1; row < src.rows - s.y; ++row ) {
            if ( src.at<uchar>(row, s.x) == obj_id ) {
                p = Point(s.x, row);
                dividingLine.push_back(p);
                searchSpace.push_back(p);
            } else {
                break;
            }
        }

        neighborhood(dividingLine, searchSpace);
        for (auto v = searchSpace.begin(); v != searchSpace.end(); ++v ) {
            p = *v;
            buffer.at<uchar>(p.y, p.x) = found;
            if ( p.x < s.x ) {
                leftSearch.push_back(p);
            } else if ( p.x > s.x ) {
                rightSearch.push_back(p);
            }
        }
    }

    vector<Point> leftVectorized = floodFillAndVectorize(
                                            src, buffer, leftSearch, found);
    vector<Point> rightVectorized = floodFillAndVectorize(
                                            src, buffer, rightSearch, found);

    // add initial point to end of left vector to maintain order
    rightVectorized.push_back(centroid(dividingLine));

    // merge vectors in order
    vector<Point> merged;
    merged.reserve(leftVectorized.size() + rightVectorized.size() );
    reverse(rightVectorized.begin(), rightVectorized.end());
    merged.insert(merged.end(), rightVectorized.begin(), rightVectorized.end());

    // added before current position
    merged.insert(merged.end(), leftVectorized.begin(), leftVectorized.end());

    return merged;
}

vector<Point> floodFillAndVectorize(Mat& src, Mat& buffer,
                                    vector<Point> searchSpace, int found) {

    vector<Point> vectorized;

    int i;
    vector<Point> neighbors;
    while ( searchSpace.size() > 0 ) {
        // calculate the center point of the new search space to be used
        // as a point in the final vector of the shape
        vectorized.push_back(centroid(searchSpace));

        // find neighborhood for the current search space
        neighbors.clear();
        neighborhood(searchSpace, neighbors);

        // find new search space based on neighbors that are a part of object
        searchSpace.clear();
        Point p;
        for ( i = 0; i < neighbors.size(); ++i ) {
            p = neighbors[i];
            if ( inBounds(p, src.cols, src.rows) ) {
                if ( src.at<uchar>(p.y, p.x) == found &&
                            buffer.at<uchar>(p.y, p.x) != found ) {
                    buffer.at<uchar>(p.y, p.x) = found;
                    searchSpace.push_back(p);
                }
            }
        }
    }

    return vectorized;
}

//http://stackoverflow.com/a/20947961
struct comparePoints {
    bool operator()(const Point & a, const Point & b) {
        return (a.x != b.x || a.y != b.y);
    }
};

void neighborhood(vector<Point> points, vector<Point>& neighbors) {
    set<Point, comparePoints> uNeighbors; // unique neighbors
    Point p;
    for ( int i = 0; i < points.size(); ++i ) {
        p = points[i];
        // above/below
        uNeighbors.insert(Point(p.x, p.y-1));
        uNeighbors.insert(Point(p.x, p.y+1));
        // right neighborhood
        uNeighbors.insert(Point(p.x-1, p.y-1));
        uNeighbors.insert(Point(p.x-1, p.y));
        uNeighbors.insert(Point(p.x-1, p.y+1));
        // left neighborhood
        uNeighbors.insert(Point(p.x+1, p.y+1));
        uNeighbors.insert(Point(p.x+1, p.y));
        uNeighbors.insert(Point(p.x+1, p.y-1));
    }
    set<Point,comparePoints>::iterator it;
    for ( it = uNeighbors.begin(); it != uNeighbors.end(); ++it ) {
        neighbors.push_back(*it);
    }
}

bool inBounds(Point p, int w, int h) {
    return (p.x >= 0 && p.x < w && p.y >= 0 && p.y < h);
}

Point borderPoint(Mat& src, int obj_id) {
    Point p;
    int row, col;
    for ( row = 0; row < src.rows; ++row ) {
        for ( col = 0; col < src.cols; ++col ) {
            if (src.at<uchar>(row, col) == obj_id) {
                p = Point(col, row);
                return p;
            }
        }
    }
    return p;
}

