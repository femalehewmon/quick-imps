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

vector<string> filesInFolder(string& dirpath) {
    dirpath.append("/*");
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

    int n = 1;                  // number of frames to stores
    vector<Mat> history;        // stores the previous n frames

    Mat orig, background, buffer, output;

    // process all video frames
    for ( int i = 0; i < img_files.size(); ++i ) {
        orig = imread(img_files[i], CV_LOAD_IMAGE_COLOR);
        if(!orig.data ) {
            cout << "No image data" <<  endl;
            return -1;
        }

        if ( !background.empty() ) {
            subtract(background, orig, output);
            cvtColor(output, output, CV_BGR2GRAY);
            threshold(output, output, 25, 255, 0);

            // show output frame
            imshow("Output", output);

            // If ESC key pressed, exit program
            if (waitKey(30) == 27) {
                break;
            }
        }

        // update the history, and merge if enough frames have passed
        history.push_back(orig.clone());
        if ( history.size() > n ) {

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
