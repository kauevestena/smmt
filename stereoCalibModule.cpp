/* This is sample from the OpenCV book. The copyright notice is below */

/* *************** License:**************************
   Oct. 3, 2008
   Right to use this code in any way you want without warranty, support or any guarantee of it working.

   BOOK: It would be nice if you cited it:
   Learning OpenCV: Computer Vision with the OpenCV Library
     by Gary Bradski and Adrian Kaehler
     Published by O'Reilly Media, October 3, 2008

   AVAILABLE AT:
     http://www.amazon.com/Learning-OpenCV-Computer-Vision-Library/dp/0596516134
     Or: http://oreilly.com/catalog/9780596516130/
     ISBN-10: 0596516134 or: ISBN-13: 978-0596516130

   OPENCV WEBSITES:
     Homepage:      http://opencv.org
     Online docs:   http://docs.opencv.org
     Q&A forum:     http://answers.opencv.org
     Issue tracker: http://code.opencv.org
     GitHub:        https://github.com/Itseez/opencv/

   ************************************************** */

/*

this code have some deep modifications by: Kaue de Moraes Vestena
 AND the GitHub OPENCV community
*/

#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/core/core.hpp"

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include <fstream>

using namespace cv;
using namespace std;

//

int ipow(int base, int exp)
{
    int result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}

//structure to decide the group of parameters for the calibration
struct calib_options
{
//bool k1k2;
//bool p1p2;
//bool k3;
//bool thin_prism;
//bool k4k5k6;
//bool tau

    int nParams,paramCombination,paramCombinationMono;

    calib_options(bool K1,bool K2,bool P1P2,bool K3,bool THIN_PRISM,bool K4K5K6,bool TAU,bool fixMonoCalib,bool fixAllMono);
};

calib_options::calib_options(bool K1,bool K2,bool P1P2,bool K3,bool THIN_PRISM,bool K4K5K6,bool TAU,bool fixMonoCalib,bool fixAllMono)
{
    nParams = 4;
    paramCombination = CALIB_USE_INTRINSIC_GUESS;

    //if(K1K2){nParams += 2;}
    //if(P1P2){nParams += 2;}
    if(K3)
    {
        nParams = 5;
    }
    if(K4K5K6)
    {
        nParams = 8;
        paramCombination +=CALIB_RATIONAL_MODEL;
    }
    if(THIN_PRISM)
    {
        nParams = 12;
        paramCombination += CALIB_THIN_PRISM_MODEL;
    }
    if(TAU)
    {
        nParams = 14;
        paramCombination += CALIB_TILTED_MODEL;
    }

    if(!K3)
    {
        paramCombination += CALIB_FIX_K3;
    }
    if(!K2)
    {
        paramCombination+=CALIB_FIX_K2;
    }
    if(!K1)
    {
        paramCombination+=CALIB_FIX_K1;
    }
    if(!P1P2)
    {
        paramCombination+=CALIB_ZERO_TANGENT_DIST;
    }

    if(!K4K5K6 && nParams > 5)
    {
        paramCombination += (CALIB_FIX_K4+CALIB_FIX_K5+CALIB_FIX_K6);
    }

    //to have the difference between the two
    paramCombinationMono = paramCombination;

    //now the only for two cameras
    if(fixMonoCalib)
    {
        paramCombination+=(CALIB_FIX_FOCAL_LENGTH+CALIB_FIX_PRINCIPAL_POINT);
    }
    if(fixAllMono)
    {
        paramCombination+=CALIB_FIX_INTRINSIC;
    }
}

// TODO (kauevestena#1#): relatory of the reprojection error per point
void computeReprojectionErrors( const vector<vector<Point3f> >& objectPoints,
                                const vector<vector<Point2f> >& imagePoints,
                                const vector<Mat>& rvecs, const vector<Mat>& tvecs,
                                const Mat& cameraMatrix, const Mat& distCoeffs,string filename)
{
    vector<Point2f> imagePoints2;
    int i,j, totalPoints = 0;
    double totalErr = 0, err,err2=0;
    //perViewErrors.resize(objectPoints.size());
    ofstream report(filename);
    Point2f temp;

    for( i = 0; i < (int)objectPoints.size(); ++i )
    {

        projectPoints( Mat(objectPoints[i]), rvecs[i], tvecs[i], cameraMatrix,
                       distCoeffs, imagePoints2);
        int n = (int)objectPoints[i].size();

        for(j=0; j<n; ++j)
        {
            //cout<<" i "<<i<<" j "<<j<<endl;

            err = cv::norm(Mat(imagePoints[i]).at<float>(i,j), Mat(imagePoints2).at<float>(i,j), CV_L2,noArray());


            report << err<<",";

            err2 += err;
        }

//        err = norm(Mat(imagePoints[i]), Mat(imagePoints2), CV_L2);

        report <<endl;

//            cout<<Mat(imagePoints[i]).at<float>(0,0)<<endl;
//            cout<<Mat(imagePoints2).at<float>(0,0)<<endl;

        //cout<<std::sqrt(err)<<endl;

        //int n = (int)objectPoints[i].size();
//        perViewErrors[i] = (float) std::sqrt(err2*err2/n);
//        //cout<<perViewErrors[i]<<endl;
//        totalErr        += err2*err2;
//        totalPoints     += n;
    }

    report.close();

//    return std::sqrt(totalErr/totalPoints);
}


string subsCharacter(string input,string character)
{
    if (input.find(character) !=-1)
    {
        return input.substr(input.find(character)+1,string::npos);
    }
    else
    {
        return input;
    }
}

static void
StereoCalib(calib_options options,const vector<string>& imagelist, Size boardSize,int npyr=0,double w=0,double h=0,const float squareSize = 1.f, bool useCalibrated=true, bool showRectified=true,bool displayCorners = true)
{
    if( imagelist.size() % 2 != 0 )
    {
        cout << "Error: the image list contains odd (non-even) number of elements\n";
        return;
    }


    ///
    //when the pyramids are aplied, we will need to apply a scale to the measurements
    int measSc = ipow(2,npyr);
    ///

    //cout<<measSc<<endl;

    //true;
    const int maxScale = 2;
    // ARRAY AND VECTOR STORAGE:

    vector<vector<Point2f> > imagePoints[2];
    vector<vector<Point3f> > objectPoints;
    Size imageSize;

    int i, j, k, nimages = (int)imagelist.size()/2;

    imagePoints[0].resize(nimages);
    imagePoints[1].resize(nimages);
    vector<string> goodImageList;

    namedWindow("corners",0);

    for( i = j = 0; i < nimages; i++ )
    {
        for( k = 0; k < 2; k++ )
        {
            const string& filename = imagelist[i*2+k];
            Mat img = imread(filename, 0);

//            /// ~-------------~ THE PYRAMID STUFF **********
//            for (int j=0; j<npyr; j++)
//            {
//                pyrDown(img,img);
//            }
//            /// ~-------------~

            if(img.empty())
                break;
            if( imageSize == Size() )
            {
                imageSize = img.size();
                cout<<Size()<<endl;
            }
            else if( img.size() != imageSize )
            {
                cout << "The image " << filename << " has the size different from the first image size. Skipping the pair\n";
                break;
            }
            bool found = false;
            vector<Point2f>& corners = imagePoints[k][j];
//            for( int scale = 1; scale <= maxScale; scale++ )
//            {
            Mat timg;
//                if( scale == 1 )
            timg = img;

            /// ~-------------~ THE PYRAMID STUFF **********
            for (int j=0; j<npyr; j++)
            {
                pyrDown(timg,timg);
            }
            /// ~-------------~
//                else
//                    resize(img, timg, Size(), scale, scale);
            found = findChessboardCorners(timg, boardSize, corners,
                                          CALIB_CB_ADAPTIVE_THRESH | CALIB_CB_NORMALIZE_IMAGE);
//                if( found )
//                {
//                    if( scale > 1 )
//                    {
//                        Mat cornersMat(corners);
//
//                        if (npyr > 0)
//                        {
//                            cornersMat *= measSc;
//                        }


//            if(!corners.empty())
//            {
//                cout<<corners.at(0)<<endl;
//            }

            if (npyr > 0)
            {

                for (unsigned int i = 0; i<corners.size(); i++)
                {
                    corners.at(i) *= float(measSc);
                }
            }

//            if(!corners.empty())
//            {
//                cout<<corners.at(0)<<endl;
//            }
//                        cout<<corners.at(0)<<endl<<endl;

//                        cornersMat *= 1./scale;
//                    }
//                    break;
//                }

//                cout<<scale<<endl;
//            }
            if( displayCorners )
            {
                cout << filename << endl;
                Mat cimg, cimg1;
                cvtColor(img, cimg, COLOR_GRAY2BGR);
                drawChessboardCorners(cimg, boardSize, corners, found);
                double sf = 640./MAX(img.rows, img.cols);
                resize(cimg, cimg1, Size(), sf, sf);
                imshow("corners", cimg1);
                char c = (char)waitKey(500);
                if( c == 27 || c == 'q' || c == 'Q' ) //Allow ESC to quit
                    exit(-1);
            }
            else
                putchar('.');
            if( !found )
                break;
            cornerSubPix(img, corners, Size(11,11), Size(-1,-1),
                         TermCriteria(TermCriteria::COUNT+TermCriteria::EPS,
                                      30, 0.01));

//                                      cout<<corners.at(corners.size()-1);
        }
        if( k == 2 )
        {
            goodImageList.push_back(imagelist[i*2]);
            goodImageList.push_back(imagelist[i*2+1]);
            j++;
        }

        //destroyAllWindows();
    }

    int nimages2 = (int)goodImageList.size()/2;

    cout << j << " pairs have been successfully detected.\n";
    nimages = j;
    if( nimages < 2 )
    {
        cout << "Error: too little pairs to run the calibration\n";
        return;
    }

    imagePoints[0].resize(nimages);
    imagePoints[1].resize(nimages);
    objectPoints.resize(nimages);


    //creating the object points
    for( i = 0; i < nimages; i++ )
    {
        for( j = 0; j < boardSize.height; j++ )
            for( k = 0; k < boardSize.width; k++ )
                objectPoints[i].push_back(Point3f(k*squareSize, j*squareSize, 0));
        //cout << objectPoints[i]<<endl;
    }

    cout << "Running stereo calibration ...\n";

    Mat cameraMatrix[2], distCoeffs[2];

/////KMs
    distCoeffs[0] = Mat::zeros(1,options.nParams,CV_32F);
    distCoeffs[1] = Mat::zeros(1,options.nParams,CV_32F);

//with that you can control what parameters will be estimated (to exclude the k3 from the game, for example)
///end

    Mat R, T, E, F;


    double rmsL,rmsR;

    ///KAUE'S MODs

    TermCriteria term(TermCriteria::COUNT+TermCriteria::EPS, 100, 1e-5);

////    cout<<endl<<imagePoints[0].at(1)<<endl<<endl;
////

//    for (unsigned int i = 0; i<imagePoints[0].size(); i++)
//    {
//        cout<<endl<<imagePoints[0].at(i)<<endl<<endl;
//
////    cout<<endl<<imagePoints[1].at(i)<<endl<<endl;
//    }


//    cout<<endl<<imagePoints[1].at(imagePoints[1].size()-1)<<endl<<endl;


    cameraMatrix[0] = initCameraMatrix2D(objectPoints,imagePoints[0],imageSize,0);
    cameraMatrix[1] = initCameraMatrix2D(objectPoints,imagePoints[1],imageSize,0);

    cout<<"Initialized Left Camera Matrix "<<cameraMatrix[0]<<endl;
    cout<<"Initialized RIght Camera Matrix "<<cameraMatrix[1]<<endl<<endl;

    vector<Mat> rotL,rotR,trL,trR;
    Mat stdInL,stdInR,stdExL,stdExR,PerViewL,PerViewR;

//    ///***
//    cout<<endl<<distCoeffs[0].size()<<endl;
//
//    cout<<endl<<distCoeffs[1].size()<<endl;
//    ///***

    //doing the calibration of both cameras separately, at first

    //the
    rmsL = calibrateCamera(objectPoints,imagePoints[0],imageSize,cameraMatrix[0],
                           distCoeffs[0],rotL,trL,stdInL,stdExL,PerViewL,options.paramCombinationMono,term);
    cout << "Single calibration of left camera done with "<<rmsL<<" RMS" <<endl;

    rmsR = calibrateCamera(objectPoints,imagePoints[1],imageSize,cameraMatrix[1],
                           distCoeffs[1],rotR,trR,stdInR,stdExR,PerViewR,options.paramCombinationMono,term);
    cout << "Single calibration of right camera done with "<<rmsR<<" RMS" <<endl;

    cout<<"Second initialized Left Camera Matrix "<<cameraMatrix[0]<<endl;
    cout<<"Second initialized RIght Camera Matrix "<<cameraMatrix[1]<<endl;

    // save intrinsic parameters of the single-camera-calibration
    FileStorage fs1("intrinsicsLeft.yml", FileStorage::WRITE);
    if( fs1.isOpened() )
    {
        fs1 << "M" << cameraMatrix[0] << "D" << distCoeffs[0];
        fs1.release();
    }

    FileStorage fs2("intrinsicsRight.yml", FileStorage::WRITE);
    if( fs2.isOpened() )
    {
        fs2 << "M" << cameraMatrix[1] << "D" << distCoeffs[1];
        fs2.release();
    }


//    vector<float> perViewErrors[2];



    computeReprojectionErrors(objectPoints,imagePoints[0],rotL,trL,cameraMatrix[0],distCoeffs[0],"reprojErrLeft.csv");
    computeReprojectionErrors(objectPoints,imagePoints[1],rotR,trR,cameraMatrix[1],distCoeffs[1],"reprojErrRight.csv");

//    cout<<" Average reprojection error in Left: "<<errL<<endl<<" Average reprojection error in Right: "<<errR<<endl;

/// END OF KAUE'S MODs

    double rms = stereoCalibrate(objectPoints, imagePoints[0], imagePoints[1],
                                 cameraMatrix[0], distCoeffs[0],
                                 cameraMatrix[1], distCoeffs[1],
                                 imageSize, R, T, E, F,
                                 options.paramCombination
                                 ,term );
    cout << "done with RMS error=" << rms << endl;

// CALIBRATION QUALITY CHECK
// because the output fundamental matrix implicitly
// includes all the output information,
// we can check the quality of calibration using the
// epipolar geometry constraint: m2^t*F*m1=0
    double err = 0;
    double pairAvgError = 0;
    double maxPairError = 0;
    double maxCjPointError = 0;
    int npt2,bestPairs=0;

    //FileStorage fs3("epipolar_report.xml", FileStorage::WRITE);
    ofstream epiData("epipolarErrors.csv");
    ofstream epiData2("epipolarErrors2.csv");
    ofstream before("distortedPoints.csv");
    ofstream after("undistortedPoints.csv");
    ofstream finalOnes("bestPairs.txt");
    ofstream report("report.txt");
    ofstream reportMono("reportExtra.txt");

    int npoints = 0;
    vector<Vec3f> lines[2];
    for( i = 0; i < nimages; i++ )
    {
        int npt = (int)imagePoints[0][i].size();
        if (i==0)
        {
            npt2=npt;
        }
        Mat imgpt[2];
        for( k = 0; k < 2; k++ )
        {
            imgpt[k] = Mat(imagePoints[k][i]);
            before << imgpt[k] <<endl;
            undistortPoints(imgpt[k], imgpt[k], cameraMatrix[k], distCoeffs[k], Mat(), cameraMatrix[k]);
            after << imgpt[k] <<endl;
            computeCorrespondEpilines(imgpt[k], k+1, F, lines[k]);
        }
        for( j = 0; j < npt; j++ )
        {
            double errij = fabs(imagePoints[0][i][j].x*lines[1][j][0] +
                                imagePoints[0][i][j].y*lines[1][j][1] + lines[1][j][2]) +
                           fabs(imagePoints[1][i][j].x*lines[0][j][0] +
                                imagePoints[1][i][j].y*lines[0][j][1] + lines[0][j][2]);

            Point3d Lpoint,Rpoint;

            Lpoint.x = imagePoints[0][i][j].x;
            Lpoint.y = imagePoints[0][i][j].y;
            Lpoint.z = 1;

            Rpoint.x = imagePoints[1][i][j].x;
            Rpoint.y = imagePoints[1][i][j].y;
            Rpoint.z = 1;


            double errij2 = sampsonDistance(Mat(Lpoint),Mat(Rpoint),F);

            ///SOME NAIVE TESTS
            ////            Mat mmm1,mmm2;
            ////            mmm1 = Mat(Lpoint);
            ////            transpose(mmm1,mmm2);
            ////
            ////            Mat multip =  mmm2* F * Mat(Rpoint);
            ////            double errij2 = multip.at<double>(0,0);
            ///END OF IT =P


            if (errij > maxCjPointError)
            {
                maxCjPointError=errij;
            }
            err += errij;
            pairAvgError += errij;
            epiData << errij <<",";
            epiData2 << errij2 <<",";
            //cout<<errij<<" numero "<<j<<endl;
            //fs3 <<"err_cj"+to_string(j+1)<< errij;
        }
        epiData<<endl;
        epiData2<<endl;
        npoints += npt;
        cout << "pair "<<npoints/npt<<" average error "<< pairAvgError/npt <<" ,pair wich left is: "<<goodImageList.at(i*2)<<" and right: "<<goodImageList.at(i*2+1)<<endl;
        if (pairAvgError/npt > maxPairError)
        {
            maxPairError=pairAvgError/npt;
        }
        //cout <<" npoints "<<npoints<<" npt "<<npt<<" i "<<i<<endl;
        if ( pairAvgError/npt < 1 )
        {
            //to generate the file with only the best stereopairs
            finalOnes << goodImageList.at(i*2)<<endl;
            finalOnes << goodImageList.at(i*2+1)<<endl;
            bestPairs++;
        }
        pairAvgError=0;
    }
    cout << "Average epipolar err = " <<  err/npoints << endl;
    cout << "Major pair average error= " <<  maxPairError << endl;
    cout << "Major single conjugate point error= " <<  maxCjPointError << endl;
    cout<< "total of best Stereopairs: " <<bestPairs<<endl;

    report << "STEREOCALIBRATION done with RMS error         " << rms << endl;
    report << "Single calibration of left camera done with the RMS         "<<rmsL<<endl;
    report << "Single calibration of right camera done with the RMS        "<<rmsR<<endl;

    report << "Average epipolar err =                                   " <<  err/npoints << endl;
    report << "Major pair average error: " <<  maxPairError << endl;
    report << "Major single conjugate point error: " <<  maxCjPointError << endl;
    report<< "total of best Stereopairs: " <<bestPairs<<endl;

    if (w != 0 && h!= 0)
    {
        //in metric units, of the final calibration values
        double fovXl,fovYl,fovXr,fovYr,Fl,Fr,Propl,Propr;
        Point2d ppL,ppR;
        //Left camera:
        calibrationMatrixValues(cameraMatrix[0],imageSize,w,h,fovXl,fovYl,Fl,ppL,Propl);
        calibrationMatrixValues(cameraMatrix[1],imageSize,w,h,fovXr,fovYr,Fr,ppR,Propr);

        report << endl;
        report <<"LEFT:"<<endl;
        report <<"Focal Distance:"<<Fl<<" mm"<<endl;
        report <<"Principal Point:"<<ppL<<" (mm) origin upper left corner"<<endl;
        report <<"FOV:"<<fovXl<<" x "<<fovYl<<" deg."<<endl;

        report << endl;
        report <<"Right:"<<endl;
        report <<"Focal Distance:"<<Fr<<" mm"<<endl;
        report <<"Principal Point:"<<ppR<<" (mm) origin upper left corner"<<endl;
        report <<"FOV:"<<fovXr<<" x "<<fovYr<<" deg."<<endl;

    }

    //The report With New Features

    reportMono<<"Câmera ESQUERDA"<<endl<<endl;

    reportMono<<"Desvios dos Intrínsecos"<<endl;
    reportMono<<stdInL<<endl;
    reportMono<<endl;

    reportMono<<"Desvios dos Extrínsecos"<<endl;
    reportMono<<stdExL<<endl;
    reportMono<<endl;

    reportMono<<"Erro Médio de Reprojeção, por frame"<<endl;
    reportMono<<PerViewL<<endl;

    reportMono<<endl<<endl<<endl;



    reportMono<<"Câmera DIREITA"<<endl<<endl;

    reportMono<<"Desvios dos Intrínsecos"<<endl;
    reportMono<<stdInR<<endl;
    reportMono<<endl;

    reportMono<<"Desvios dos Extrínsecos"<<endl;
    reportMono<<stdExR<<endl;
    reportMono<<endl;

    reportMono<<"Erro Médio de Reprojeção, por frame"<<endl;
    reportMono<<PerViewR<<endl;


    reportMono<<endl;


    epiData.close();
    before.close();
    after.close();
    finalOnes.close();
    report.close();
    epiData2.close();
    reportMono.close();

    //fs3.release();

    //end of calibration stuff, now its time to rectification

    // save intrinsic parameters
    FileStorage fs("intrinsics.yml", FileStorage::WRITE);
    if( fs.isOpened() )
    {
        fs << "M1" << cameraMatrix[0] << "D1" << distCoeffs[0] <<
           "M2" << cameraMatrix[1] << "D2" << distCoeffs[1];
        fs.release();
    }
    else
        cout << "Error: can not save the intrinsic parameters\n";

    Mat R1, R2, P1, P2, Q;
    Rect validRoi[2];

    stereoRectify(cameraMatrix[0], distCoeffs[0],
                  cameraMatrix[1], distCoeffs[1],
                  imageSize, R, T, R1, R2, P1, P2, Q,
                  CALIB_ZERO_DISPARITY, 1, imageSize, &validRoi[0], &validRoi[1]);

    fs.open("relative_extrinsics.yml", FileStorage::WRITE);
    if( fs.isOpened() )
    {
        fs <<"F"<<F<<"E"<<E<< "R" << R << "T" << T << "R1" << R1 << "R2" << R2 << "P1" << P1 << "P2" << P2 << "Q" << Q;
        fs.release();
    }
    else
        cout << "Error: can not save the extrinsic parameters\n";

    // OpenCV can handle left-right
    // or up-down camera arrangements
    //bool isVerticalStereo = fabs(P2.at<double>(1, 3)) > fabs(P2.at<double>(0, 3));
    bool isVerticalStereo = false;
    //cout<<isVerticalStereo<<endl;

// COMPUTE AND DISPLAY RECTIFICATION
    //Precompute maps for cv::remap()

    Mat rmap[2][2];

    initUndistortRectifyMap(cameraMatrix[0], distCoeffs[0], R1, P1, imageSize, CV_16SC2, rmap[0][0], rmap[0][1]);
    initUndistortRectifyMap(cameraMatrix[1], distCoeffs[1], R2, P2, imageSize, CV_16SC2, rmap[1][0], rmap[1][1]);
    Mat canvas,imgT;
    double sf;

    if( !isVerticalStereo )
    {
        sf = 600./MAX(imageSize.width, imageSize.height);
        w = cvRound(imageSize.width*sf);
        h = cvRound(imageSize.height*sf);
        canvas.create(h, w*2, CV_8UC3);
    }
    else
    {
        sf = 300./MAX(imageSize.width, imageSize.height);
        w = cvRound(imageSize.width*sf);
        h = cvRound(imageSize.height*sf);
        canvas.create(h*2, w, CV_8UC3);
    }

    for( i = 0; i < nimages; i++ )
    {
        for( k = 0; k < 2; k++ )
        {
            Mat img = imread(goodImageList[i*2+k], 0), rimg, cimg;
//            cout<<subsCharacter(goodImageList[i*2+k],"/")<<endl;
//            cout<<k<<endl;

///KAUE's MODs

            undistort(img,imgT,cameraMatrix[k],distCoeffs[k]);
            if(npyr < 2)
            {pyrDown(imgT,imgT);
            pyrDown(imgT,imgT);}
            namedWindow("undistorted",0);
            imshow("undistorted",imgT);
            imwrite("undistorted/u"+subsCharacter(goodImageList[i*2+k],"/"),imgT);
            imgT.release();

///END OF IT

            remap(img, rimg, rmap[k][0], rmap[k][1], INTER_LINEAR);
            cvtColor(rimg, cimg, COLOR_GRAY2BGR);
            Mat canvasPart = !isVerticalStereo ? canvas( Rect(w*k, 0, w, h)) : canvas(Rect(0, h*k, w, h));
            resize(cimg, canvasPart, canvasPart.size(), 0, 0, INTER_AREA);
            if( useCalibrated )
            {
                Rect vroi(cvRound(validRoi[k].x*sf), cvRound(validRoi[k].y*sf),
                          cvRound(validRoi[k].width*sf), cvRound(validRoi[k].height*sf));
                rectangle(canvasPart, vroi, Scalar(0,0,255), 3, 8);
            }
        }

        //drawing the epilines
        if( !isVerticalStereo )
            for( j = 0; j < canvas.rows; j += 16 )
                line(canvas, Point(0, j), Point(canvas.cols, j), Scalar(0, 255, 0), 1, 8);
        else
            for( j = 0; j < canvas.cols; j += 16 )
                line(canvas, Point(j, 0), Point(j, canvas.rows), Scalar(0, 255, 0), 1, 8);
        imshow("rectified", canvas);
        char c = (char)waitKey();
        if( c == 27 || c == 'q' || c == 'Q' )
            break;
    }
}


static bool readStringList( const string& filename, vector<string>& l )
{
    l.resize(0);
    FileStorage fs(filename, FileStorage::READ);
    if( !fs.isOpened() )
        return false;
    FileNode n = fs.getFirstTopLevelNode();
    if( n.type() != FileNode::SEQ )
        return false;
    FileNodeIterator it = n.begin(), it_end = n.end();
    for( ; it != it_end; ++it )
        l.push_back((string)*it);
    return true;
}

//vector::string stereoImgList(string pathLeft,string pathRight)
//{
//    vector::string left,right,res;
//
//
//    return res;
//}

int main()
{
//samples:
//    Size boardSize = Size(9,6);
//    string imagelistfn = "list2.xml";

//my data:
    string imagelistfn = "list.xml";
    string pathL = "/home/kauevestena/data/calib31102016/esquerda/*.JPG";
    string pathR = "/home/kauevestena/data/calib31102016/direita/*.JPG";

    //originalRES

    vector<string> imagelist;
    bool ok = readStringList(imagelistfn, imagelist);
    if(!ok || imagelist.empty())
    {
        cout<<"creating the imagelist file"<<endl;

        cout<<"run again to use the file"<<endl;

        ofstream imgList("list.xml");

        //header
        imgList<<"<?xml version=\"1.0\"?>"<<endl;
        imgList<<"<opencv_storage>"<<endl;
        imgList<<"<images>"<<endl;

        //list of images

        std::vector<cv::String> listL,listR;

        cv::glob(pathL,listL,false);
        cv::glob(pathR,listR,false);

        if (listL.size() != listR.size())
        {
            cout<<"há um número diferente de imagens na pasta da esquerda e da direita"<<endl;
            cout<<"considerando a esquerda, há  "<<listL.size()-listR.size()<<" imagens a mais nela"<<endl;

            return -1;

        }
        else if (listL.size() == 0 || listR.size() == 0)
        {
            cout<<"uma das duas pastas se encontra vazia..."<<endl;

            return -1;
        }
        else
        {
            cout<<"recording the path to the images"<<endl<<endl;

            for (unsigned int i = 0; i<listL.size(); i++)
            {
                imgList<<listL.at(i)<<endl;
                imgList<<listR.at(i)<<endl;

                cout<<i<<",";
            }

            cout<<endl<<"total de "<<listL.size()<<" estereopares"<<endl;
        }

        //footer
        imgList<<"</images>"<<endl;
        imgList<<"</opencv_storage>"<<endl;

        imgList.close();

        return 0;
    }

    cout<<endl;

//    for (unsigned int i = 0;i<imagelist.size();i++)
//    {
//        cout<<imagelist.at(i)<<endl;
//    }

//cout << subsCharacter("teste/testao","/")<<endl;

    calib_options the_options(1,1,1,0,0,0,0,0,1);
    //bool K1,bool K2,bool P1P2,bool K3,bool THIN_PRISM,bool K4K5K6,bool TAU,bool fixMonoCalib,bool fixAllMono


    //novo tabuleiro
    Size boardSize = Size(6,5);
    //velho tabuleiro
//    Size boardSize = Size(7,8);

    StereoCalib(the_options,imagelist, boardSize,2,22.3,14.9,7.8921);


    return 0;
}
