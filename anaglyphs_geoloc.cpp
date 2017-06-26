#include <stdio.h>
#include <iostream>
#include <stdio.h>
#include <iostream>

#include <string>
#include <vector>

#include "opencv2/core.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/xfeatures2d.hpp"
#include "opencv2/calib3d.hpp"
#include "opencv2/imgproc.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <typeinfo>
#include <vector>
#include <ostream>
#include <iomanip>


using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;

void purgeEntries(String textToSearch,vector<String> *input)
{
//  vector<String> temp = *input;
    int purged = 0;
    for (unsigned int i = 0; i < input->size(); i++)

//  vector<string>::iterator it;
//
//  for (it = input->begin();it != input->end;it++)
    {
        if (input->at(i).find(textToSearch) == -1)
        {
            input ->erase(input->begin()+i);
            purged++;
        }
    }

//    for (unsigned int i = 0;i<temp.size();i++)
//    {
//        if (temp.at(i).find(textToSearch) == -1)
//        {
//            temp.erase()
//        }
//    }

//    cout<<"purged "<<purged<<" entries"<<endl;
//    *input = temp;
}

int numberFromString(int nPlaces,String input)
{
    String substring = input.substr(input.size()-4-nPlaces,nPlaces);

//    cout<<substring<<endl;

    return atoi(substring.c_str());
}

vector<string> SepString(string input,char defSeparator,vector<string> separators)
{
    vector<string> res;

//    cout<<input<<endl;

//    cout<<input<<endl;

    int pos = 0;

    //replace all non-default with the first one
    while (pos != -1)
    {
        for (unsigned int i = 0; i < separators.size(); i++)
        {
            pos = input.find(separators.at(i));
//            cout<<separators.at(i)<<",";

            if(pos > -1)
            {
                input.replace(pos,1,&defSeparator);
            }

        }
    }
//    cout<<endl;

    //now split off the separators
    std::istringstream split(input);
    //all elegant: https://goo.gl/TpgLBM --- thank you user 560648
    for (std::string temp; std::getline(split, temp, defSeparator); res.push_back(temp));

//    if (res.size()>23 && res.at(0) == "#MARK1PVAA" )
//    {cout<<input<<endl;nops = true;}

    return res;
}

vector<vector<string>> logreaderBasic(string logname)
{
    vector<vector<string>> res;

    vector<string> temp;

    std::ifstream data;
    data.open(logname);
    string line;

    //any non-comma field separator must be HERE
    vector<string> separators = {" "};
    //comma is the default separator, but if it was changed, change here too
    char separator = ',';


//    unsigned line_count = 1;

    while (std::getline(data,line))
    {

//        string::iterator end_it = utf8::find_invalid(line.begin(), line.end());
//        if (end_it != line.end()) {
//            cout << "Invalid UTF-8 encoding detected at line " << line_count << "\n";
//            cout << "This part is fine: " << string(line.begin(), end_it) << "\n";
//            line_count++;
//        }

        temp = SepString(line,separator,separators);
//        if (nops)
//        {cout<<line;nops=false;}

        res.push_back(temp);

        temp.clear();
    }

    data.close();

    return res;
}

Mat filterMatches(const std::vector<cv::DMatch> matches,const std::vector<cv::KeyPoint>&keypoints1,const std::vector<cv::KeyPoint>& keypoints2,std::vector<cv::DMatch>& goodMatches,double dist,double confidence_interval,int op)
{
    goodMatches.clear();
    // converting keypoints to just 2D points
    std::vector<cv::Point2f> points1, points2;
    for (std::vector<cv::DMatch>::const_iterator it= matches.begin(); it!= matches.end(); ++it)
    {
        // from de matches, extract just the keypoints of the left frame
        float x= keypoints1[it->queryIdx].pt.x;
        float y= keypoints1[it->queryIdx].pt.y;
        points1.push_back(cv::Point2f(x,y));

        // from de matches, extract just the keypoints of the left frame
        x= keypoints2[it->trainIdx].pt.x;
        y= keypoints2[it->trainIdx].pt.y;
        points2.push_back(cv::Point2f(x,y));
    }
    // calculating the Fundamental mat with the ransac algorithm
    std::vector<uchar> inliers(points1.size(),0);

    Mat res;

    if (op == 1)
    {
        res = findFundamentalMat(cv::Mat(points1),cv::Mat(points2),inliers,CV_FM_RANSAC,dist,confidence_interval);
    }
    else
    {
        res = findHomography(cv::Mat(points1),cv::Mat(points2),inliers,RANSAC,dist);
    }

//        cv::Mat fundamental=         Mat H =
    // with the inliers, choose from the matches only the "good" ones
    std::vector<uchar>::const_iterator
    itIn= inliers.begin();
    std::vector<cv::DMatch>::const_iterator
    itM= matches.begin();
    // for all matches
    for ( ; itIn!= inliers.end(); ++itIn, ++itM)
    {
        if (*itIn)
        {
            // it is a valid match
            goodMatches.push_back(*itM);
        }
    }

    return res;
}

void keypoints2points(std::vector<cv::Point2f> &points1,std::vector<cv::Point2f> &points2,const std::vector<cv::KeyPoint>&keypoints1,const std::vector<cv::KeyPoint>&keypoints2,const std::vector<cv::DMatch> matches)
{
    for (std::vector<cv::DMatch>::const_iterator it= matches.begin(); it!= matches.end(); ++it)
    {
        float x,y;

        // from the matches, extract just the keypoints of the left frame
        x= keypoints1[it->queryIdx].pt.x;
        y= keypoints1[it->queryIdx].pt.y;
        points1.push_back(cv::Point2f(x,y));

        // from the matches, extract just the keypoints of the right frame
        x= keypoints2[it->trainIdx].pt.x;
        y= keypoints2[it->trainIdx].pt.y;
        points2.push_back(cv::Point2f(x,y));
    }
}

void anaglyphFromPair(String leftPath,String rightPath,string destPath = "out.jpg",bool show = false)
{
    Mat left = imread(leftPath);

    Mat right = imread(rightPath);
//    Mat teste2 = teste;

    Ptr<Feature2D> features = ORB::create();

    vector<KeyPoint> kpointsl,kpointsr;
    Mat descripl,descripr;

    features -> detect(left,kpointsl);

    features -> detect(right,kpointsr);

    features -> compute(left,kpointsl,descripl);

    features -> compute(right,kpointsr,descripr);

//    cout<<kpoints.size();

//    drawKeypoints(left,kpointsl,left,Scalar::all(-1),1);

    vector<DMatch> matches,goodMatches,goodMatches2;

    BFMatcher matcher(NORM_HAMMING,true);

    matcher.match(descripl,descripr,matches);

    Mat imMatches;

//    cout<<kpointsl.size()<<endl;cout<<kpointsr.size()<<endl;

    Mat F,H;

    //two-step filtering

    F = filterMatches(matches,kpointsl,kpointsr,goodMatches,3,.99,1);

    H = filterMatches(goodMatches,kpointsl,kpointsr,goodMatches2,3,.99,2);

    //transforming the right image
    Mat rTrans;// = right;
//
//    double tx = H.at<double>(0,2);
//    double ty = H.at<double>(1,2);
//    double h[3][3] = {{1, 0, tx}, {0, 1, ty}, {0, 0, 1}};
//    H = Mat(3, 3, CV_64F, h);

//    cout<<H<<endl;

    warpPerspective(right,rTrans,H,right.size(),INTER_LINEAR+WARP_INVERSE_MAP);

    Mat anag;// = rTrans;


    //splitting images
    Mat splitLeft[3];
    Mat splitRight[3];

    split(left,splitLeft);
    split(rTrans,splitRight);

    vector<Mat> channels;


    channels.push_back(splitRight[1]);
    channels.push_back(splitRight[2]);
    channels.push_back(splitLeft[0]);

    cv::merge(channels,anag);

//    cout<<kpointsl.size()<<endl;cout<<kpointsr.size()<<endl;

//    cout<<matches.size()<<endl<<goodMatches.size()<<endl<<goodMatches2.size()<<endl;
//
//    cout<<F<<endl<<H<<endl;

//    cout << <<endl<< H.at<double>(1,2)<<endl;

    if (show)
    {
        drawMatches(left,kpointsl,right,kpointsr,goodMatches,imMatches);

        namedWindow("teste",0);
        namedWindow("teste2",0);

        imshow("teste",imMatches);

        imshow("teste2",anag);

        waitKey();
    }

    imwrite(destPath,anag);


}

struct stereoPairData
{
    String leftPath,rightPath;

    double lat,lgt,h,yaw,pitch,roll;

    int idx;

    void printData();
};

void stereoPairData::printData()
{
    cout<<lat<<","<<lgt<<","<<h<<","<<yaw<<","<<pitch<<","<<roll<<endl;
}

struct allStereopairData
{
    vector<stereoPairData> allData;

    vector<string> anaglyphNames;

    void createAnaglyphs(void);

    void createKML(void);

    allStereopairData(string pathToPoseData,string pathToLPhotos,string pathToRPhotos);
};

void allStereopairData::createAnaglyphs()
{
    for (unsigned int i = 0; i < allData.size(); i++)
    {
        string anaglyphName = "imgs/anag"+to_string(allData.at(i).idx)+".jpg";
        cout<<anaglyphName<<endl;
        //anaglyphFromPair(allData.at(i).leftPath,allData.at(i).rightPath,anaglyphName);
        anaglyphNames.push_back(anaglyphName);
    }

}

void allStereopairData::createKML()
{
    system("mkdir imgs");

    ofstream out("anaglifos.kml");

    //thank you http://articles.extension.org/sites/default/files/w/7/70/GPS_Benchmark.kmz
    //for the model in which this file have been based

    out<<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>"<<endl;
    out<<"<kml xmlns=\"http://earth.google.com/kml/2.1\">"<<endl;
    out<<"<Document>"<<endl;
    out<<"<name>anaglifos_SMMT</name>"<<endl;
    out<<"<LookAt>"<<endl;
    out<<"<longitude>"<<allData.at(0).lgt<<"</longitude>"<<endl;
    out<<"<latitude>"<<allData.at(0).lat<<"</latitude>"<<endl;
    out<<"<altitude>"<<allData.at(0).h<<"</altitude>"<<endl;
    out<<"<range>100.00</range>"<<endl;
    out<<"<tilt>0.00</tilt>"<<endl;
    out<<"<heading>0.00</heading>"<<endl;
    out<<"</LookAt>"<<endl;
    out<<"<StyleMap id=\"000000_0_0_copy66\">"<<endl;
    out<<"<Pair>"<<endl;
    out<<"<key>normal</key>"<<endl;
    out<<"<styleUrl>#000000_0_0_copy67</styleUrl>"<<endl;
    out<<"</Pair>"<<endl;
    out<<"<Pair>"<<endl;
    out<<"<key>highlight</key>"<<endl;
    out<<"<styleUrl>#000000_0_0_copy68</styleUrl>"<<endl;
    out<<"</Pair>"<<endl;
    out<<"</StyleMap>"<<endl;
    out<<"<Style id=\"000000_0_0_copy68\">"<<endl;
    out<<"<IconStyle>"<<endl;
    out<<"<scale>0.7</scale>"<<endl;
    out<<"<Icon>"<<endl;
    out<<"<href>http://maps.google.com/mapfiles/kml/paddle/ylw-diamond.png</href>"<<endl;
    out<<"</Icon>"<<endl;
    out<<"<hotSpot x=\"32\" y=\"1\" xunits=\"pixels\" yunits=\"pixels\"/>"<<endl;
    out<<"</IconStyle>"<<endl;
    out<<"<LabelStyle>"<<endl;
    out<<"<color>ff00ffff</color>"<<endl;
    out<<"</LabelStyle>"<<endl;
    out<<"<ListStyle>"<<endl;
    out<<"<ItemIcon>"<<endl;
    out<<"<href>http://maps.google.com/mapfiles/kml/paddle/ylw-diamond-lv.png</href>"<<endl;
    out<<"</ItemIcon>"<<endl;
    out<<"</ListStyle>"<<endl;
    out<<"</Style>"<<endl;
    out<<"<Style id=\"000000_0_0_copy67\">"<<endl;
    out<<"<IconStyle>"<<endl;
    out<<"<scale>0.7</scale>"<<endl;
    out<<"<Icon>"<<endl;
    out<<"<href>http://maps.google.com/mapfiles/kml/paddle/ylw-diamond.png</href>"<<endl;
    out<<"</Icon>"<<endl;
    out<<"<hotSpot x=\"32\" y=\"1\" xunits=\"pixels\" yunits=\"pixels\"/>"<<endl;
    out<<"</IconStyle>"<<endl;
    out<<"<LabelStyle>"<<endl;
    out<<"<color>ff00ffff</color>"<<endl;
    out<<"</LabelStyle>"<<endl;
    out<<"<ListStyle>"<<endl;
    out<<"<ItemIcon>"<<endl;
    out<<"<href>http://maps.google.com/mapfiles/kml/paddle/ylw-diamond-lv.png</href>"<<endl;
    out<<"</ItemIcon>"<<endl;
    out<<"</ListStyle>"<<endl;
    out<<"</Style>"<<endl;
    out<<"<Style id=\"000000_0_0_copy74\">"<<endl;
    out<<"<IconStyle>"<<endl;
    out<<"<scale>0.7</scale>"<<endl;
    out<<"<Icon>"<<endl;
    out<<"<href>http://maps.google.com/mapfiles/kml/paddle/ylw-diamond.png</href>"<<endl;
    out<<"</Icon>"<<endl;
    out<<"<hotSpot x=\"32\" y=\"1\" xunits=\"pixels\" yunits=\"pixels\"/>"<<endl;
    out<<"</IconStyle>"<<endl;
    out<<"<LabelStyle>"<<endl;
    out<<"<color>ff00ffff</color>"<<endl;
    out<<"</LabelStyle>"<<endl;
    out<<"<ListStyle>"<<endl;
    out<<"<ItemIcon>"<<endl;
    out<<"<href>http://maps.google.com/mapfiles/kml/paddle/ylw-diamond-lv.png</href>"<<endl;
    out<<"</ItemIcon>"<<endl;
    out<<"</ListStyle>"<<endl;
    out<<"</Style>"<<endl;
    out<<"<Style id=\"000000_0_0_copy78\">"<<endl;
    out<<"<IconStyle>"<<endl;
    out<<"<scale>0.7</scale>"<<endl;
    out<<"<Icon>"<<endl;
    out<<"<href>http://maps.google.com/mapfiles/kml/paddle/ylw-diamond.png</href>"<<endl;
    out<<"</Icon>"<<endl;
    out<<"<hotSpot x=\"32\" y=\"1\" xunits=\"pixels\" yunits=\"pixels\"/>"<<endl;
    out<<"</IconStyle>"<<endl;
    out<<"<LabelStyle>"<<endl;
    out<<"<color>ff00ffff</color>"<<endl;
    out<<"</LabelStyle>"<<endl;
    out<<"<ListStyle>"<<endl;
    out<<"<ItemIcon>"<<endl;
    out<<"<href>http://maps.google.com/mapfiles/kml/paddle/ylw-diamond-lv.png</href>"<<endl;
    out<<"</ItemIcon>"<<endl;
    out<<"</ListStyle>"<<endl;
    out<<"</Style>"<<endl;
    out<<"<Style id=\"000000_0_0_copy63\">"<<endl;
    out<<"<IconStyle>"<<endl;
    out<<"<scale>0.7</scale>"<<endl;
    out<<"<Icon>"<<endl;
    out<<"<href>http://maps.google.com/mapfiles/kml/paddle/ylw-diamond.png</href>"<<endl;
    out<<"</Icon>"<<endl;
    out<<"<hotSpot x=\"32\" y=\"1\" xunits=\"pixels\" yunits=\"pixels\"/>"<<endl;
    out<<"</IconStyle>"<<endl;
    out<<"<LabelStyle>"<<endl;
    out<<"<color>ff00ffff</color>"<<endl;
    out<<"</LabelStyle>"<<endl;
    out<<"<ListStyle>"<<endl;
    out<<"<ItemIcon>"<<endl;
    out<<"<href>http://maps.google.com/mapfiles/kml/paddle/ylw-diamond-lv.png</href>"<<endl;
    out<<"</ItemIcon>"<<endl;
    out<<"</ListStyle>"<<endl;
    out<<"</Style>"<<endl<<endl;

    out.precision(10);

    //placemark creation
    for (unsigned int i = 0; i < allData.size(); i++)
    {
        out  <<  "<Placemark>" <<endl;
//        out  <<  "<name maxLines=\"1\">" <<  to_string(i+1)  << "</name>" <<endl;
out  <<  "<name maxLines=\"1\">" <<  " " << "</name>" <<endl;
        out  <<  "<description><![CDATA[<table width=\"358\" border=\"0\">" <<endl;
        out  <<  "<tr>" <<endl;
        out  <<  "<td>"  << "gerado com dados do SMMT do LAPE, acesse o " <<  "<a href=\"https://lapeufpr.blogspot.com.br/\">blog do LAPE</a> </td>" <<endl;
        out  <<  "</tr>" <<endl;
        out  <<  "<tr>" <<endl;
        out  <<  "<td><IMG SRC=\""<<  anaglyphNames.at(i)  <<"\" WIDTH=\""<<  5184/8  <<"\" HEIGHT=\""<<  3456/8  <<"\" BORDER=\"1\" ALT=\"platform\"></td>" <<endl;
        out  <<  "</tr>" <<endl;
        out  <<  "</table>]]></description>" <<endl;
        out  <<  "<LookAt>" <<endl;
        out  <<  "<longitude>" <<  allData.at(i).lat  << "</longitude>" <<endl;
        out  <<  "<latitude>" <<  allData.at(i).lgt  << "</latitude>" <<endl;
        out  <<  "<altitude>" <<  allData.at(i).h  << "</altitude>" <<endl;
        out  <<  "<range>100.00</range>" <<endl;
        out  <<  "<tilt>0.00</tilt>" <<endl;
        out  <<  "<heading>0.00</heading>" <<endl;
        out  <<  "</LookAt>" <<endl;
        out  <<  "<styleUrl>#000000_0_0_copy78</styleUrl>" <<endl;
        out  <<  "<Point>" <<endl;
        out  <<  "<coordinates>"<<  allData.at(i).lgt  <<","<<  allData.at(i).lat <<","<<  allData.at(i).h  <<"</coordinates>" <<endl;
        out  <<  "</Point>" <<endl;
        out  <<  "</Placemark>" <<endl<<endl;
    }

    out<<"</Document>"<<endl<<"</kml>";

}


allStereopairData::allStereopairData(string pathToPoseData,string pathToLPhotos,string pathToRPhotos)
{
    vector<vector<string>> poses = logreaderBasic("poses.txt");

    vector<String> origNamesL,origNamesR;

    cv::glob(pathToLPhotos,origNamesL,false);
    cv::glob(pathToRPhotos,origNamesR,false);

    //purge all non-jpg files
    purgeEntries(".jpg",&origNamesL);
    purgeEntries(".jpg",&origNamesR);

    //    for (unsigned int i )

    //go through all images
    int major_list  = numberFromString(6,origNamesL.at(origNamesL.size()-1));
    int major_list2 = numberFromString(6,origNamesR.at(origNamesR.size()-1));

    if (major_list2 < major_list)
    {
        major_list = major_list2;

    }

    int lastValid = -1;
    for (unsigned int i = 0; i < major_list; i++)
    {
        int indexLeft  = numberFromString(6,origNamesL.at(i));
        int indexRight = numberFromString(6,origNamesR.at(i));



        if (indexLeft == indexRight && indexLeft == i + 1)
        {


            stereoPairData temp;
            temp.leftPath  = origNamesL.at(i);
            temp.rightPath = origNamesR.at(i);
            temp.lat   = strtod(poses.at(i).at(0).c_str(),nullptr);
            temp.lgt   = strtod(poses.at(i).at(1).c_str(),nullptr);
            temp.h     = strtod(poses.at(i).at(2).c_str(),nullptr);
            temp.pitch = strtod(poses.at(i).at(3).c_str(),nullptr);
            temp.roll  = strtod(poses.at(i).at(4).c_str(),nullptr);
            temp.yaw   = strtod(poses.at(i).at(5).c_str(),nullptr);
            temp.idx = i;

//            temp.printData();

            allData.push_back(temp);
            lastValid = indexLeft;
        }
        else
        {
            break;
            //TODO: implement the logic to handle the absence of a valid pair
        }

    }

    createAnaglyphs();
             createKML();


}


int main( int argc, char** argv )
{
//             string leftPath  = "/home/kauevestena/data/undistorted/esquerda/left000155.jpg";
//             string rightPath = "/home/kauevestena/data/undistorted/direita/right000155.jpg";
//
//             anaglyphFromPair(leftPath,rightPath);

    //the position and orientations forall photos
    //    vector<vector<string>> poses = logreaderBasic("poses.txt");

    allStereopairData job("poses.txt",
                          "/home/kauevestena/data/undistorted/direita/","/home/kauevestena/data/undistorted/esquerda/");


    return 0;
}



