#include <opencv2/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <string>
#include <vector>
#include <opencv2/core/persistence.hpp>
#include "exif.h"
#include "exif.cpp"

using namespace cv;

double exif_read(std::string img_path,bool print=false)
{

    double seconds;

    double h,m,s,ms;

    // Read the JPEG file into a buffer
    FILE *fp = fopen(img_path.c_str(), "rb");
    if (!fp)
    {
        printf("Can't open file.\n");
    }

    fseek(fp, 0, SEEK_END);
    unsigned long fsize = ftell(fp);
    rewind(fp);
    unsigned char *buf = new unsigned char[fsize];
    if (fread(buf, 1, fsize, fp) != fsize)
    {
        printf("Can't read file.\n");
        delete[] buf;
    }
    fclose(fp);

    // Parse EXIF
    easyexif::EXIFInfo result;
    int code = result.parseFrom(buf, fsize);
    delete[] buf;
    if (code)
    {
        printf("Error parsing EXIF: code %d\n", code);
        return 0.0;
    }

        h = double(atoi(result.DateTimeOriginal.substr(11,2).c_str()));
        m = double(atoi(result.DateTimeOriginal.substr(14,2).c_str()));
        s = double(atoi(result.DateTimeOriginal.substr(17,2).c_str()));

        ms = double(atoi(result.SubSecTimeOriginal.c_str()));

        seconds = h*3600 + m*60 + s + ms/100;

    if (print)
    {
        cout<<result.DateTimeOriginal<<"     "<<result.SubSecTimeOriginal<<"    "<<seconds<<endl;
    }

    return seconds;
}

string preffix(unsigned int number,int nZeros)
{
    int maxNumb = 10;

    string res = "";

    for(unsigned int i = 1; i<nZeros; i++)
    {
        maxNumb *= 10;
    }

    if (maxNumb < number)
    {
        return "";
    }

    int cont = maxNumb;

    for(;;)
    {
        if (cont > number)
        {
            res+="0";
            cont/= 10;
        }
        else
        {
            break;
        }
    }

//cout<<maxNumb;

    return res;
}

//void imagelist(string path,vector<String>* output,bool outputxt = false)
//{
//
//    vector<String> fn;
//    cv::glob(path,fn,false);
//
//    ofstream file;
//    if(outputxt)
//    {
//        file.open("list.txt");
//    }
//
//    for (int i = 0; i < fn.size(); i++)
//    {
//
//        if(outputxt)
//        {
//            file << fn.at(i)<<endl;
//        }
//
//    }
//
//    *output = fn;
//}

struct photoReader
{
    //path for the original images
    string originalPath,calibPath,outPath;
    string LfolderName,RfolderName,Ext;

    bool saveImgs;

    //the calibration data
    //camera matrices and distortion arrays
    cv::Mat McamL,McamR,VdistL,VdistR;

    //vectors to store the names of the original images
    vector<String> origNamesL,origNamesR;

    //vectors to store the names of the undistorted images
    vector<string> outNamesL,outNamesR;

    //to store the EXIF timestamps
    vector<double> exifTimestampsL,exifTimestampsR,exifTimeDiffs;

    photoReader(string origPathName,string ext,string lFname,string rFname,string outFname,string calibDataPath,bool writeImgs=true);

    void readingCalibdata(bool print=false);

    void imagelist2(bool print = false);

    void outputImgNames(bool print=false);

    void rectifyImages();
};

void photoReader::readingCalibdata(bool print)
{
    FileStorage fs(calibPath, FileStorage::READ);

    fs["M1"] >> McamL;
    fs["D1"] >> VdistL;
    fs["M2"] >> McamR;
    fs["D2"] >> VdistR;

    //the scaling of the matrices should be removed before exporting the intrinsics
    //TO BE REMOVED:
    //McamL *=4;McamR*=4;

    ///debugging
    if  (print)
    {
        cout<<endl<<endl<<McamL<<endl<<endl<<VdistL<<endl<<endl<<McamR<<endl<<endl<<VdistR<<endl<<endl<<endl;
    }

    fs.release();
}

void photoReader::imagelist2(bool print)
{
    string lsep1="";
    string rsep1="";
    string sep2="";

    if(originalPath.back() != '/')
    {
        lsep1+="/";
        rsep1+="/";
    }

    if(LfolderName.back() != '/' || RfolderName.back() != '/' )
    {
        sep2+="/";
    }

    if(Ext.find(".")==-1)
    {
        sep2+=".";
        if(Ext.find("*")==-1)
        {
            sep2+="*";
        }
    }


    //LEFT CAMERA
    string toReadL = originalPath+lsep1+LfolderName+sep2+Ext;
    //RIGHT CAMERA
    string toReadR = originalPath+rsep1+RfolderName+sep2+Ext;

    if(toReadL.find("*")==-1 || toReadR.find("*")==-1)
    {
        cout<<"O programa nao sera capaz de ler  as imagens, confira os argumentos passados"<<endl;
    }

    if (print)
    {
        cout<<toReadL<<endl<<toReadR<<endl;
    }

    cv::glob(toReadL,origNamesL,false);
    cv::glob(toReadR,origNamesR,false);
    //cout<<endl<<origNamesL.size();
}

void photoReader::outputImgNames(bool print)
{

cout<<origNamesL.size()<<"  "<<origNamesR.size()<<endl;

    if (origNamesL.size() == origNamesR.size())
    {
        if (print)
        {
            for (unsigned int i = 0; i<origNamesL.size(); i++)
            {
                cout<<origNamesL.at(i)<<endl<<origNamesR.at(i)<<endl<<endl;
            }
        }
    }
    else
    {
        cerr<<"há mais imagens de uma câmera que da outra, conferir";
    }
}

void photoReader::rectifyImages()
{
    //creating the two directories
    string left  = outPath+"/esquerda";
    string right = outPath+"/direita";

    string outname;

    string cmdL = "mkdir "+left;
    string cmdR = "mkdir "+right;

    // #linux_specific
    system(cmdL.c_str());
    system(cmdR.c_str());

    cv::Mat orig,mod;

    //reading, rectifying, saving and storing the full path
    //LEFT

    //The "fakes" are a workaround for some photos that aren't taken
    for (unsigned int i = 0; i<origNamesL.size(); i++)
    {
        string fakepart = "";

        bool fakefound = false;

        //read
        if(saveImgs){
        orig = cv::imread(origNamesL.at(i),IMREAD_COLOR);}

//        cout<<origNamesL.at(i).find("fake")<<endl;

//                cout<<origNamesR.at(i)<<endl;

        if (origNamesL.at(i).find("fake") < origNamesL.at(i).size() || origNamesL.at(i).find("fake") < 0)
        {
//            cout<<origNamesL.at(i).find("fake")<<endl;
            fakepart = "_fake";
            fakefound = true;
        }

        //compensate for lens distortion
        if(saveImgs){
        undistort(orig,mod,McamL,VdistL);}

        //the outname
        outname = left+"/left"+preffix(i+1,5)+to_string(i+1)+fakepart+".jpg";
        outNamesL.push_back(outname);
        if(saveImgs){
        cout<<"undistorting and saving "<<outname<<endl;}
//        else
//        {cout<<"reading "<<outname<<endl;}

        //saving the image
        if(saveImgs){
        imwrite(outname,mod);
        orig.release();}

        outname.clear();
        mod.release();

        //the EXIF infos
        if (!fakefound)
        {
            double tempSeconds = exif_read(origNamesL.at(i));
            exifTimestampsL.push_back(tempSeconds);
        }

    }

    cout<<endl;

    //RIGHT
    for (unsigned int i = 0; i<origNamesR.size(); i++)
    {
        string fakepart = "";

        bool fakefound = false;

        //read
        if(saveImgs){
        orig = cv::imread(origNamesR.at(i),IMREAD_COLOR);}

        if (origNamesR.at(i).find("fake") < origNamesR.at(i).size() || origNamesR.at(i).find("fake") < 0)
        {

            fakepart = "_fake";
            fakefound = true;
        }

        //compensate for lens distortion
        if(saveImgs){
        undistort(orig,mod,McamR,VdistR);}
//        else
//        {cout<<"reading "<<outname<<endl;}

        //the outname
        outname = right+"/right"+preffix(i+1,5)+to_string(i+1)+fakepart+".jpg";
        outNamesR.push_back(outname);
        if(saveImgs){
        cout<<"undistorting and saving "<<outname<<endl;}


        //saving the image
        if(saveImgs){
        imwrite(outname,mod);
        orig.release();}

        outname.clear();
        mod.release();

        //the EXIF infos
        if (!fakefound)
        {
        double tempSeconds = exif_read(origNamesR.at(i));
        exifTimestampsR.push_back(tempSeconds);
        }
    }

}

photoReader::photoReader(string origPathName,string ext,string lFname,string rFname,string outFname,string calibDataPath,bool writeImgs)
{
    ///CONSTRUCTOR
    originalPath = origPathName;
    calibPath = calibDataPath;
    LfolderName = lFname;
    RfolderName = rFname;
    Ext = ext;
    outPath = outFname;
    saveImgs = writeImgs;

    //to read calibration data
    readingCalibdata();
    //to create the lists of images
    imagelist2();
    //a little test
    outputImgNames();

    //making the directory for the
    string command01 = "mkdir "+outPath;
    // #linux_specific
    system(command01.c_str());

    rectifyImages();

/// Timekeeping functions to implement
//////////    for (unsigned int i = 0;i<exifTimestampsL.size();i++)
//////////    {
//////////        int dif1 = fabs(exifTimestampsL.at(i)-exifTimestampsR.at(i));
//////////
//////////        if (dif1 > 0.1)
//////////        {
//////////            cout<<"diferença na timestamp do EXIF"<<endl;
//////////        }
//////////
//////////        if (i > 0)
//////////        {
//////////            double dif21 = fabs(exifTimestampsL.at(i) - exifTimestampsL.at(i-1));
//////////            exifTimeDiffs.push_back(dif21);
//////////        }
//////////    }

    //    cout<<exifTimeDiffs.size()<<"   "<<exifTimestampsL.size()<<endl;

}
