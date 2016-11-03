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

vector<string> filesInFolder(string& dirpath, string extra="");

vector<string> filesInFolder(string& dirpath, string extra) {
    dirpath.append(extra);
    glob_t glob_result;
    glob( dirpath.c_str(), GLOB_TILDE, NULL, &glob_result );
    vector<string> files;
    for ( int i = 0; i < glob_result.gl_pathc; ++i ) {
        files.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return files;
}

void combineHistory(Mat& curr, vector<Mat> history, Mat& dst) {
    Mat diff;
    for ( int i = 0; i < history.size(); ++i ) {
        subtract(history[i], curr, diff);
        threshold(diff, diff, 30, 255, 0);
        diff.copyTo(dst, diff);
    }
}

Mat mergedHistoryMask(vector<string> img_files) {
    int N = img_files.size() / 2; // number of frames to store
    cout << N << endl;
    vector<Mat> history;          // stores the previous N frames

    Mat morph_element = getStructuringElement(MORPH_ELLIPSE, Size(5, 5));
    Mat orig;
    for ( int i = 0; i < N; ++i ) {
        orig = imread(img_files[i], CV_LOAD_IMAGE_GRAYSCALE);
        if(!orig.data ) {
            cout << "No image data" <<  endl;
            return orig;
        }
        history.push_back(orig.clone());
    }

    Mat output = Mat::zeros(history[0].size(), history[0].type());
    Mat buffer = Mat::zeros(history[0].size(), history[0].type());
    for ( int i = 0; i < img_files.size(); ++i ) {
        orig = imread(img_files[i], CV_LOAD_IMAGE_GRAYSCALE);
        if(!orig.data ) {
            cout << "No image data" <<  endl;
            return orig;
        }

        combineHistory(orig, history, buffer);
        morphologyEx(buffer, buffer, MORPH_CLOSE, morph_element);

        Mat labels, stats, centroids;
        int n = connectedComponentsWithStats(buffer, labels, stats, centroids);
        for ( int j = 1; j < n; ++j ) {
            if ( stats.at<int>(j, CC_STAT_AREA) > 100 ) {
                // Isolate the labeled region
                compare(labels, j, buffer, CMP_EQ);
                buffer.copyTo(output, buffer);
            }
        }
    }

    return output;
}

int main(int argc, char** argv)
{
    if ( argc < 2 ) {
        cout << "ERROR: Missing required data directory argument" << endl;
        return -1;
    }

    string dirpath = argv[1];
    vector<string> img_dirs = filesInFolder(dirpath, "/*");
    for ( int i = 0; i < img_dirs.size(); ++i ) {
        vector<string> img_files = filesInFolder(img_dirs[i], "/*masked.bmp");
        if ( img_files.size() == 0 ) {
            cout << "Data files in " << dirpath << " were not found!" <<  endl;
            return -1;
        }

        Mat mergedMask = mergedHistoryMask(img_files);
        Mat orig = imread(img_files[i], CV_LOAD_IMAGE_GRAYSCALE);
        Mat display = Mat::zeros(orig.size(), orig.type());
        addWeighted(mergedMask, 0.5, orig, 0.75, 0, display);
        putText(display, img_dirs[i], Point(20, 20), 0, 0.5, 255);

        imshow("Output", display);
        waitKey(0);
    }

    return 0;
}
