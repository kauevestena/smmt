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
//#include "utf8.h"
//
//#include <locale>
//#include <codecvt>

using namespace std;

//bool nops = false;


vector<string> SepString(string input,char defSeparator,vector<string> separators)
{
    vector<string> res;

//    cout<<input<<endl;

    //replace all non-default with the first one
    for (unsigned int i = 0; i < separators.size(); i++)
    {
        size_t pos = input.find(separators.at(i));
//        cout<<pos<<",";

        if(pos > -1)
        {
            input.replace(pos,1,&defSeparator);
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

vector<vector<string>> logreaderBasic(string logname,vector<string> separators,char defSeparator)
{
    vector<vector<string>> res;

    vector<string> temp;

    std::ifstream data;
    data.open(logname);
    string line;



//	std::locale utf8_locale(locale(), new codecvt_utf8<wchar_t>);
//	data.imbue(utf8_locale);

    unsigned line_count = 1;

    while (std::getline(data,line))
    {

//        string::iterator end_it = utf8::find_invalid(line.begin(), line.end());
//        if (end_it != line.end()) {
//            cout << "Invalid UTF-8 encoding detected at line " << line_count << "\n";
//            cout << "This part is fine: " << string(line.begin(), end_it) << "\n";
//            line_count++;
//        }

        temp = SepString(line,defSeparator,separators);
//        if (nops)
//        {cout<<line;nops=false;}

        res.push_back(temp);

        temp.clear();
    }

    data.close();

//cout<<res.size()<<endl;

    return res;
}




struct gpsTimeB
{
    int week;
    double seconds;
};

//struct that contains the INS/GNSS observations
struct obs2
{

    gpsTimeB obsTime;

    //position
    double lat;
    double lgt;
    double h;

    //velocity
    double velE;
    double velN;
    double velU;

    //attitude
    double pitch;    //along Y axis ("phi")
    double roll;     //along X axis ("omega")
    double azimuth;  //along Z axis ("kappa")

    //status
    string status,type;
};

//struct that contains the covariance matrices of the INS/GNSS observations
struct obsCov2
{
    gpsTimeB covTime;

    double pos_mvc[3][3];
    double vel_mvc[3][3];
    double att_mvc[3][3];

};

//struct that contains all the others of the INS nested
struct obsWcov2
{
    obs2 observation;

    obsCov2 obsCovariances;
};

struct survData2
{
    bool bugstatus=false;

    std::vector<obsWcov2> pvaWcov;

    vector<obs2> pva;

    vector<obsCov2> covs;

    vector<vector<string>> data;

    //any non-comma field separator must be HERE
    vector<string> otherSeparators = {";","*"};
    //comma is the default separator, but if it was changed, change here too
    char separator = ',';

    vector<double> diffTimes;

    survData2(string logname,bool gen_log);

    void dataSplitter(void);

    void dataInterpolation(void);
};

void survData2::dataSplitter()
{
//ofstream pvaTimes("obsTimes.txt");
//ofstream covTimes("covsTimes.txt");
//ofstream outNav("navData.txt");

//because of some bug nb = no-bug
//    int nb = 0;

    size_t np = 4;

    //failsafe
    if (data.empty())
    {
        return;
    }

    for (unsigned int i = 0; i < data.size(); i++)
    {
        //one for each message
        if (data.at(i).at(0) == "#MARK1PVAA")
        {
            obs2 tempObs;

            tempObs.lat = strtod(data.at(i).at(12).c_str(),nullptr); //12


            if (fabs(tempObs.lat) > 90)
            {
                cout<<data.at(i).size()<<",";
                cout<<"bug"<<endl;
                bugstatus = true;
//                nb = 1;
            }

            tempObs.lgt = strtod(data.at(i).at(13).c_str(),nullptr); //13
            tempObs.h   = strtod(data.at(i).at(14).c_str(),nullptr); //14

            tempObs.velE = strtod(data.at(i).at(16).c_str(),nullptr); //16
            tempObs.velN = strtod(data.at(i).at(15).c_str(),nullptr); //15
            tempObs.velU = strtod(data.at(i).at(17).c_str(),nullptr); //17

            tempObs.roll    = strtod(data.at(i).at(18).c_str(),nullptr); //18
            tempObs.pitch   = strtod(data.at(i).at(19).c_str(),nullptr); //19
            tempObs.azimuth = strtod(data.at(i).at(20).c_str(),nullptr); //20

            tempObs.status = data.at(i).at(21);  //21

            tempObs.obsTime.week = atoi(data.at(i).at(10).c_str());    //10
            tempObs.obsTime.seconds = strtod(data.at(i).at(11).c_str(),nullptr);  //11

            tempObs.type = "MARKPVA";

//            pvaTimes<< setprecision(15)  << tempObs.obsTime.seconds<<"  "<<i<<endl;

//            outNav<<setprecision(15)<< tempObs.lat<<"   "<<tempObs.lgt<<"   "<<tempObs.h<<" "<<tempObs.pitch<<" "<<tempObs.roll<<"  "<<tempObs.azimuth<<endl;

            pva.push_back(tempObs);
        }

        if (data.at(i).at(0) == "%INSCOVSA")
        {
            obsCov2 tempCov;

            //attitude covariances
            tempCov.att_mvc[0][0] = strtod(data.at(i).at(5).c_str(),nullptr);
            tempCov.att_mvc[0][1] = strtod(data.at(i).at(6).c_str(),nullptr);
            tempCov.att_mvc[0][2] = strtod(data.at(i).at(7).c_str(),nullptr);

            tempCov.att_mvc[1][0] = strtod(data.at(i).at(8).c_str(),nullptr);
            tempCov.att_mvc[1][1] = strtod(data.at(i).at(9).c_str(),nullptr);
            tempCov.att_mvc[1][2] = strtod(data.at(i).at(10).c_str(),nullptr);

            tempCov.att_mvc[2][0] = strtod(data.at(i).at(11).c_str(),nullptr);
            tempCov.att_mvc[2][1] = strtod(data.at(i).at(12).c_str(),nullptr);
            tempCov.att_mvc[2][2] = strtod(data.at(i).at(13).c_str(),nullptr);

            //position covariances
            tempCov.pos_mvc[0][0] = strtod(data.at(i).at(14).c_str(),nullptr); //xx
            tempCov.pos_mvc[0][1] = strtod(data.at(i).at(15).c_str(),nullptr);
            tempCov.pos_mvc[0][2] = strtod(data.at(i).at(16).c_str(),nullptr);

            tempCov.pos_mvc[1][0] = strtod(data.at(i).at(17).c_str(),nullptr);
            tempCov.pos_mvc[1][1] = strtod(data.at(i).at(18).c_str(),nullptr); //yy
            tempCov.pos_mvc[1][2] = strtod(data.at(i).at(19).c_str(),nullptr);

            tempCov.pos_mvc[2][0] = strtod(data.at(i).at(20).c_str(),nullptr);
            tempCov.pos_mvc[2][1] = strtod(data.at(i).at(21).c_str(),nullptr);
            tempCov.pos_mvc[2][2] = strtod(data.at(i).at(22).c_str(),nullptr); //zz

            //velocity covariances
            tempCov.vel_mvc[0][0] = strtod(data.at(i).at(23).c_str(),nullptr);
            tempCov.vel_mvc[0][1] = strtod(data.at(i).at(24).c_str(),nullptr);
            tempCov.vel_mvc[0][2] = strtod(data.at(i).at(25).c_str(),nullptr);

            tempCov.vel_mvc[1][0] = strtod(data.at(i).at(26).c_str(),nullptr);
            tempCov.vel_mvc[1][1] = strtod(data.at(i).at(27).c_str(),nullptr);
            tempCov.vel_mvc[1][2] = strtod(data.at(i).at(28).c_str(),nullptr);

            tempCov.vel_mvc[2][0] = strtod(data.at(i).at(29).c_str(),nullptr);
            tempCov.vel_mvc[2][1] = strtod(data.at(i).at(30).c_str(),nullptr);
            tempCov.vel_mvc[2][2] = strtod(data.at(i).at(31).c_str(),nullptr);

            tempCov.covTime.week    = atoi(data.at(i).at(1).c_str());
            tempCov.covTime.seconds = strtod(data.at(i).at(4).c_str(),nullptr);

//            covTimes<< setprecision(15) << tempCov.covTime.seconds<<"  "<<i<<endl;

            covs.push_back(tempCov);
        }
    }

}

void survData2::dataInterpolation()
{


    unsigned int pMinN,pMinP;

    double delta;

    obsWcov2 temp;

    for(unsigned int i = 0; i<pva.size(); i++)
    {
        double minN = -1000;
        double minP = 1000;

        temp.observation = pva.at(i);

        //calculating the difference in times
        for (unsigned int j = 0; j<covs.size(); j++)
        {
            delta = pva.at(i).obsTime.seconds - covs.at(j).covTime.seconds;

//            if (pva.at(i).obsTime.seconds == 0.00 || covs.at(j).covTime.seconds == 0.00)
//            {
//                cout<<"zero detectado!!"<<endl;
//            }

            if (delta <= 0 && delta > minN)
            {
                minN = delta;
                pMinN=j;
            }

            if (delta >= 0 && delta < minP)
            {
                minP = delta;
                pMinP=j;
            }

        }
//        cout<<minN<<"   "<<minP<<"  ";

        //assingning the nearest neighbour
//        cout<<minN<<"   "<<minP<<"  "<<pMinN<<" "<<pMinP<<endl;
        if (fabs(minN) <= fabs(minP))
        {
            temp.obsCovariances = covs.at(pMinN);
//            cout<<"neg"<<endl;
        }
        else
        {
            temp.obsCovariances = covs.at(pMinP);
//            cout<<"pos"<<endl;
        }

        pvaWcov.push_back(temp);
    }
}

survData2::survData2(string logname,bool gen_log)
{
    data = logreaderBasic(logname,otherSeparators,separator);

    dataSplitter();

    dataInterpolation();

    ofstream poses;

    for (unsigned int i = 0; i<pvaWcov.size(); i++)
    {
        if (i > 0)
        {
            double tempDiff = fabs(pvaWcov.at(i).observation.obsTime.seconds - pvaWcov.at(i-1).observation.obsTime.seconds);
            diffTimes.push_back(tempDiff);
        }

        if(gen_log)
        {
//        cout<<pvaWcov.at(i).observation.h<<endl;

            if (i == 0)
            {
            poses.open("poses.txt");

            poses<<"lat,lgt,h,roll,pitch,azimuth,varPos,varOri,weektime"<<endl;
            }

            poses<<std::setprecision(15)<<pvaWcov.at(i).observation.lat<<","<<pvaWcov.at(i).observation.lgt<<","<<pvaWcov.at(i).observation.h<<",";
            poses<<pvaWcov.at(i).observation.roll<<","<<pvaWcov.at(i).observation.pitch<<","<<pvaWcov.at(i).observation.azimuth<<",";
            poses<<(pvaWcov.at(i).obsCovariances.pos_mvc[0][0]+pvaWcov.at(i).obsCovariances.pos_mvc[1][1]+pvaWcov.at(i).obsCovariances.pos_mvc[2][2])/3<<",";
            poses<<(pvaWcov.at(i).obsCovariances.att_mvc[0][0]+pvaWcov.at(i).obsCovariances.att_mvc[1][1]+pvaWcov.at(i).obsCovariances.att_mvc[2][2])/3<<",";
            poses<<pvaWcov.at(i).observation.obsTime.seconds<<endl;

        }
    }

    cout<<"total: "<<pvaWcov.size()<<"  "<<pva.size()<<endl;
    if(bugstatus)
    {
        cout<<"bug 01 encontrado!"<<endl;
    }
}
