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

vector<string> filesInFolder(string& dirpath, vector<string>& files, string extra="");

vector<string> filesInFolder(string& dirpath, vector<string>& files, string extra) {
    dirpath.append(extra);
    glob_t glob_result;
    glob( dirpath.c_str(), GLOB_TILDE, NULL, &glob_result );
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
        threshold(diff, diff, 75, 255, 0);
        diff.copyTo(dst, diff);
    }
}

Mat mergedHistoryMask(vector<string> img_files) {
    int N = img_files.size(); // number of frames to store
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
    vector<string> img_dirs;
    filesInFolder(dirpath, img_dirs, "/*");
    vector<string> img_files;
    for ( int i = 0; i < img_dirs.size(); ++i ) {
        filesInFolder(img_dirs[i], img_files, "/*masked.bmp");
        if ( img_files.size() == 0 ) {
            cout << "Data files in " << dirpath << " were not found!" <<  endl;
            return -1;
        }

        cout << img_dirs[i] << endl;

        Mat mergedMask = mergedHistoryMask(img_files);
        Mat output = Mat::zeros(mergedMask.size(), mergedMask.type());

        Mat labels, stats, centroids;
        int n = connectedComponentsWithStats(mergedMask, labels, stats, centroids);
        for ( int j = 1; j < n; ++j ) {
            if ( stats.at<int>(j, CC_STAT_AREA) > 400 ) {
                // Isolate the labeled region
                compare(labels, j, mergedMask, CMP_EQ);
                mergedMask.copyTo(output, mergedMask);
            }
        }

        // first attempt: look for hole in motion mask
        Mat inverted;
        bitwise_not(output, inverted);
        output = Scalar(0);
        n = connectedComponentsWithStats(inverted, labels, stats, centroids);
        for ( int j = 1; j < n; ++j ) {
            if ( stats.at<int>(j, CC_STAT_AREA) < 5000 &&
                  stats.at<int>(j, CC_STAT_AREA) > 50 ) {
                cout << stats.at<int>(j, CC_STAT_AREA) << endl;
                // Isolate the labeled region
                compare(labels, j, inverted, CMP_EQ);
                inverted.copyTo(output, inverted);
            }
        }

        // second attempt: if smallish motion region found at the top,
        // take the bottom of that region down to first found lower region
        // third attempt: look for left or right of large full area region

        Mat orig = imread(img_files[0], CV_LOAD_IMAGE_GRAYSCALE);
        Mat display = Mat::zeros(orig.size(), orig.type());
        addWeighted(output, 0.5, orig, 0.75, 0, display);

        putText(display, img_dirs[i], Point(20, 20), 0, 0.5, 255);

        imshow("Output", display);
        waitKey(0);
        img_files.clear();
    }

    return 0;
}
