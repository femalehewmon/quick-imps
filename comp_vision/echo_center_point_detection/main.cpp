#include <cv.h>
#include <highgui.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/video/background_segm.hpp"
#include "opencv2/video/video.hpp"
#include "opencv2/video/tracking.hpp"
#include <opencv2/core/core.hpp>

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <glob.h>

using namespace cv;
using namespace std;

vector<string> filesInFolder(string& dirpath);
void mergeHistory(vector<Mat> history, Mat& dst);
bool inBounds(Rect r1, Rect r2, int dist=0);
Rect mergeRects(Rect r1, Rect r2);

vector<string> filesInFolder(string& dirpath) {
    dirpath.append("/*masked.bmp");
    glob_t glob_result;
    glob( dirpath.c_str(), GLOB_TILDE, NULL, &glob_result );
    vector<string> files;
    for ( int i = 0; i < glob_result.gl_pathc; ++i ) {
        files.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return files;
}

void mergeHistory(vector<Mat> history, Mat& dst) {
    int rows = history[0].rows;
    int cols = history[0].cols;
    int channels = history[0].channels();

    vector<uchar> merging;
    int i, row, col;
    for ( row = 0; row < rows; ++row ) {
        for ( col = 0; col < cols*channels; ++col ) {
            merging.clear();
            for ( i = 0; i < history.size(); ++i ) {
                merging.push_back(history[i].at<uchar>(row, col));
            }
            // take the median intensity across all history frames for
            // this single pixel
            size_t n = merging.size() / 2;
            nth_element(merging.begin(), merging.begin()+n, merging.end());
            int median = merging[n];
            dst.at<uchar>(row, col) = merging[n];

            // alternatively, take the mean
            //double avg = mean(merging)[0];
            //dst.at<uchar>(row, col) = avg;
        }
    }
}

void combineHistory(Mat& curr, vector<Mat> history, Mat& dst) {
    dst = Mat::zeros(curr.size(), CV_8UC1);

    Mat diff;
    for ( int i = 0; i < history.size(); ++i ) {
        subtract(history[i], curr, diff);
        cvtColor(diff, diff, CV_BGR2GRAY);
        threshold(diff, diff, 35, 255, 0);
        diff.copyTo(dst, diff);
    }
}

bool inBounds(Rect r1, Rect r2, int dist) {
    return ((r2.x < (r1.x + r1.width + dist) && (r2.x + r2.width) > (r1.x - dist)) &&
            (r2.y < (r1.y + r1.height + dist) && (r2.y + r2.height) > (r1.y - dist)) );
}

Rect mergeRects(Rect r1, Rect r2) {
    int x = r1.x < r2.x ? r1.x : r2.x;
    int y = r1.y < r2.y ? r1.y : r2.y;
    int w = r1.x + r1.width > r2.x + r2.width ? (r1.x + r1.width - x) : (r2.x + r2.width - x);
    int h = r1.y + r1.height > r2.y + r2.height ? (r1.y + r1.height - y) : (r2.y + r2.height - y);
    Rect merged = Rect(x, y, w, h);
    return merged;
}

struct RectCompare {
    bool operator() (const Rect& r1, const Rect& r2) const {
        return (r1.x < r2.x && r1.y < r2.y);
    }
};

int main(int argc, char** argv)
{
    if ( argc < 2 ) {
        cout << "ERROR: Missing required data directory argument" << endl;
        return -1;
    }

    string dirpath = argv[1];
    vector<string> img_files = filesInFolder(dirpath);
    if ( img_files.size() == 0 ) {
        cout << "Data files in " << dirpath << " were not found!" <<  endl;
        return -1;
    }

    int N = 10;                  // number of frames to store
    vector<Mat> history;         // stores the previous N frames

    Mat morph_element = getStructuringElement(MORPH_ELLIPSE, Size(5, 5));
    Mat orig, background, buffer, output, display;

    // process all video frames
    for ( int i = 0; i < img_files.size(); ++i ) {
        orig = imread(img_files[i], CV_LOAD_IMAGE_COLOR);
        if(!orig.data ) {
            cout << "No image data" <<  endl;
            return -1;
        }
        display = Mat::zeros(orig.size(), orig.type());
        buffer = Mat::zeros(orig.size(), orig.type());

        if ( !background.empty() ) {
            //subtract(background, orig, buffer);
            //cvtColor(buffer, buffer, CV_BGR2GRAY);
            //threshold(buffer, buffer, 25, 255, 0);
            combineHistory(orig, history, buffer);

            // try to merge nearby areas with morphology
            morphologyEx(buffer, buffer, MORPH_CLOSE, morph_element);

            Mat labels, stats, centroids;
            int n = connectedComponentsWithStats(
                    buffer, labels, stats, centroids);

            int close_to_label = 0;
            map<Rect, vector<int>, RectCompare> close_to_connected;
            for ( int j = 1; j < n; ++j ) {
                if ( stats.at<int>(j, CC_STAT_AREA) > 100 ) {
                    // Isolate the labeled region
                    compare(labels, j, buffer, CMP_EQ);
                    buffer.copyTo(output, buffer);

                    // Create bounding box of labeled area
                    Rect rect = boundingRect(buffer);

                    vector<int> connected;
                    // See if bounding box falls within bounds of existing
                    // bounding box. If it does, merge the regions. Otherwise,
                    // add new bounding box to list.
                    /*
                    vector<Rect> merge;
                    for ( auto const& box : close_to_connected ) {
                        if ( inBounds(box.first, rect) ) {
                            merge.push_back(box.first);
                        }
                    }

                    while ( merge.size() > 0 ) {
                        if ( merge.size() > 0 ) {
                            for ( int k = 0; k < merge.size(); ++k ) {
                                // add overlapping box labels to new entry
                                connected.insert(connected.end(),
                                             close_to_connected[merge[k]].begin(),
                                             close_to_connected[merge[k]].end());
                                // remove overlapping box from final list
                                close_to_connected.erase(merge[k]);
                                // merge overlapping boxes and add to new entry
                                rect = mergeRects(rect, merge[k]);
                            }
                        }

                        merge.clear();
                        for ( auto const& box : close_to_connected ) {
                            if ( inBounds(box.first, rect) ) {
                                merge.push_back(box.first);
                            }
                        }
                    }

                    */
                    // if no overlapping boxes found, create new list
                    if ( connected.size() == 0 ) {
                        connected.push_back(j);
                    }

                    // add final rect to list of final bounding boxes
                    close_to_connected[rect] = connected;
                }
            }

            // Draw final bounding boxes
            /*
            for ( auto const& box : close_to_connected ) {
                Rect rect = box.first;
                Point p1(rect.x, rect.y);
                Point p2(rect.x+rect.width, rect.y+rect.height);
                rectangle(output, p1, p2, 255, 1);
            }
            */

            // show output frame
            cvtColor(output, buffer, CV_GRAY2BGR);
            //addWeighted(orig, 0.5, buffer, 1.0, 0, display);
            addWeighted(buffer, 0.5, orig, 0.75, 0, display);

            imshow("Output", display);
            //output = Scalar(0);
            //display = Scalar(0);
            //buffer = Scalar(0);

            // if ESC key pressed, exit program
            if (waitKey(30) == 27) {
                break;
            }
            // if end of video reached, loop
            i = (i == img_files.size() - 1) ? 0 : i;
        }

        // update the history, and merge if enough frames have passed
        history.push_back(orig.clone());
        if ( history.size() > N ) {

            // cycle history by deleting oldest entry
            history.erase(history.begin());

            // merge history into a single Mat
            background = Mat::zeros(orig.size(), orig.type());
            mergeHistory(history, background);
        }
    }

    waitKey(0);
    return 0;
}
