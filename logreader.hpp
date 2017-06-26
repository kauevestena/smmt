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
// to the time conversions (do not forget the files!!):
#include "time_conversion.c"
#include "time_conversion.h"

/***

    THIS IS THE DEPRECATED VERSION, PLEASE, USE THE FUNCTIONS AND STRUCTS OF THE FILE
                                "logreader2.hpp"

***/

//#include <boost/filesystem.hpp>


///
const int toUTC = 3;
///

using namespace std;

//const long double pi = 4*atan(1);
//const int gpsLibLeap = 14;


struct gpsTime
{
	unsigned short week;
	double seconds;
};

//alternative for the standart definition
struct gpsTime2
{
	int week;
	double seconds;

};

//struct that contains the INS/GNSS observations
struct obs
{

	gpsTime obsTime;

	//position
	double lat;
	double lgt;
	double h;

	//velocity
	double velX;
	double velY;
	double velZ;

	//attitude
	double pitch;    //along Y axis ("phi")
	double roll;     //along X axis ("omega")
	double azimuth;  //along Z axis ("kappa")

	//status
	string status;
};

//struct that contains the covariance matrices of the INS/GNSS observations
struct obsCov
{
	gpsTime covTime;

	double pos_mvc[3][3];
	double vel_mvc[3][3];
	double att_mvc[3][3];

};

//struct that contains all the others of the INS nested
struct obsWcov
{
	obs observation;

	obsCov obsCovariances;
};

//struct for the phototimestamp (actually of the stereopair)

struct moment
{
double hour;
double minute;
double second;
};

struct stereopairTimestamp
{
double tsyear;
double tsmonth;
double tsday;

double leapSeconds;

moment leftStart;
moment leftEnd;
moment RightStart;
moment RightEnd;

gpsTime _leftStart;
gpsTime _leftEnd;
gpsTime _RightStart;
gpsTime _RightEnd;

gpsTime finalTime;
};

//struct for the closest observations to the timestamps of the photos
struct intervals
{
double shortestNeg; int shortestNegPos;
double shortestPos; int shortestPosPos;
double theTime;
};

struct survData
{    //store the INS/GNSS data in:
	std::vector<obsWcov> pvaWcov;
	//the iterator:
	//std::vector<obsWcov>::iterator  pvaWcovIt;


	//store the stereopair timestamps
	std::vector<stereopairTimestamp> timestamps;
	//std::vector<stereopairTimestamp>::iterator timestampsIt;

	//store the shortests intervals of time
	std::vector<intervals> finalIntervals;

	    //store the final ones:
	std::vector<obsWcov> finalObsWcovs;

	string log;

	string log2_;
	double year;
	double month;
	double day;
	double LeapSeconds;

	survData(string NavLog,string PhotoTimesLog,double Year,double Month,double Day,double LeapSeconds_);

	void logReader();

	void photoTimestampReader();

	void dataInterpolation();

	};

	survData::survData(string NavLog,string PhotoTimesLog,double Year,double Month,double Day,double LeapSeconds_)
	{
        log = NavLog;
        log2_ = PhotoTimesLog;
        year = Year;
        month = Month;
        day = Day;
        LeapSeconds = LeapSeconds_;

        logReader();
        photoTimestampReader();
        dataInterpolation();

	}

//
//    //store the INS/GNSS data in:
//	std::vector<obsWcov> pvaWcov;
//	//the iterator:
//	//std::vector<obsWcov>::iterator  pvaWcovIt;
//
//
//	//store the stereopair timestamps
//	std::vector<stereopairTimestamp> timestamps;
//	//std::vector<stereopairTimestamp>::iterator timestampsIt;
//
//	//store the shortests intervals of time
//	std::vector<intervals> finalIntervals;
//
//	    //store the final ones:
//	std::vector<obsWcov> finalObsWcovs;

double linTerp(double time,double timeBefore,double timeAfter,double valBefore,double valAfter)
{
//the classic linear interpolation
double value = valBefore + (  (valAfter-valBefore)* ((time-timeBefore)/(timeAfter-timeBefore)) );

return value;
}

void survData::logReader()
{
    //3 steps, to pull of the separators
	std::ifstream infile2(log);
	string line2;
    ofstream arquivo;
	arquivo.open("parte1.txt");
	while (std::getline(infile2, line2,','))
	{
    arquivo<<line2<<endl;
	}
    infile2.close();
    arquivo.close();

    ifstream parte2;
    parte2.open("parte1.txt");
    ofstream parte2b("parte2.txt");
    string line3;
    //parte2b.open("parte2.txt");
    while (std::getline(parte2,line3,';'))
	{parte2b << line3<<endl;}
    parte2b.close(); parte2.close();
    //a.close();

    ifstream parte3;
    parte3.open("parte2.txt");
    ofstream parte3b;
    string line4;
    parte3b.open("parte3.txt");
    while (std::getline(parte3,line4,'*'))
	{parte3b << line4<<endl;}
	parte3.close();
	parte3b.close();

	//now, its possible to do the data split simpler than before
	ifstream data; data.open("parte3.txt");
	string value;
	bool itsObs=false,itsCovs=false;
	bool almostOne = false; //no observation without covariance
	int cObs=0,cCovs=0;
	obs obTemp;
	obsCov covTemp;
	obsWcov bothTemp;

	ofstream report;report.open("report.txt");


        while (std::getline(data,value))
        {
            //for the observations:
            if (value=="%INSPVASA" && almostOne)
            { //loks for an observation
                itsObs = true;
                cObs++;
                continue;
            }
            if(cObs==1 && itsObs)
            {
                obTemp.obsTime.week=stod(value);
                //obTemp.obsTime.week=stringToDouble(value);
                // //cout<<"semana: "<<stod(value)<<"   "<<value<<"   "<<endl;
                cObs++;continue;
            }
            if(cObs==2 && itsObs)
            {
                cObs++;continue;
            }
            if(cObs==3 && itsObs)
            {
                cObs++;continue;
            }
            if(cObs==4 && itsObs)
            {
                obTemp.obsTime.seconds=stod(value);
                //cout<<"segundos: "<<stringTdouble2(value)<<"   "<<value<<"   "<<endl;
                //report<<value.find('.');
                //report<<stringTdouble2(value);
                //report << StringToNumber<double> ( value ) <<endl;
                cObs++;continue;
            }
            if(cObs==5 && itsObs)
            {
                obTemp.lat=stod(value);
                //cout<<"latitude: "<<value<<endl;
                cObs++;continue;
            }
            if(cObs==6 && itsObs)
            {
                obTemp.lgt=stod(value);
                //cout<<"longitude: "<<value<<endl;
                cObs++;continue;
            }
            if(cObs==7 && itsObs)
            {
                obTemp.h=stod(value);
               // cout<<"h: "<<value<<endl;
                cObs++;continue;
            }
            if(cObs==8 && itsObs)
            {
                obTemp.velX=stod(value);
                //cout<<"velX: "<<value<<endl;
                cObs++;continue;
            }
                        if(cObs==9 && itsObs)
            {
                obTemp.velY=stod(value);
                //cout<<"velY: "<<value<<endl;
                cObs++;continue;
            }
                        if(cObs==10 && itsObs)
            {
                obTemp.velZ=stod(value);
                //cout<<"velZ: "<<value<<endl;
                cObs++;continue;
            }
                    if(cObs==11 && itsObs)
            {
                obTemp.roll=stod(value);
                //cout<<"roll: "<<value<<endl;
                cObs++;continue;
            }
                if(cObs==12 && itsObs)
            {
                obTemp.pitch=stod(value);
                //cout<<"pitch: "<<value<<endl;
                cObs++;continue;
            }
            if(cObs==13 && itsObs)
            {
                obTemp.azimuth=stod(value);
                //cout<<"azimuth: "<<value<<endl;
                cObs++;continue;
            }
            if(cObs==14 && itsObs)
            {
                obTemp.status=value;
                //cout<<"status "<<value<<endl;
                cObs++;continue;
            }
            if(cObs==15 && itsObs)
            {
           // cout<<endl<<endl<<"FIM DE OBSERVAÇÃO"<<endl;
            bothTemp.observation=obTemp;
            bothTemp.obsCovariances=covTemp;
            pvaWcov.push_back(bothTemp);
            report <<std::setprecision(15)<< bothTemp.observation.obsTime.seconds<<"  "<<bothTemp.obsCovariances.covTime.seconds<<endl;
            //pva.push_back(obTemp);
            itsObs=false;cObs=0;continue;
            }


            //for the covariances
            if (value=="%INSCOVSA")
            { //loks for a covariance
                if (!almostOne) {almostOne=true;}
                itsCovs = true;
                cCovs++;
                continue;
            }

            if(cCovs==1 && itsCovs)
            {
                covTemp.covTime.week=stod(value);
                //covTemp.covTime.week=stringToDouble(value);
           ////     cout<<"semana: "<<stod(value)<<"   "<<value<<endl;
                cCovs++;continue;
            }
            if(cCovs==2 && itsCovs)
            {
            ////    cout<<"segundos: "<<stod(value)<<"   "<<value<<endl;
                covTemp.covTime.seconds=stod(value);
                cCovs++;continue;
            }
            if(cCovs==3 && itsCovs)
            {
                cCovs++;continue;
            }
            if(cCovs==4 && itsCovs)
            {
                cCovs++;continue;
            }
            if(cCovs==5 && itsCovs)
            {
                covTemp.pos_mvc[0][0]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }
                        if(cCovs==6 && itsCovs)
            {
                covTemp.pos_mvc[0][1]=stod(value);
           ////     cout<<value<<endl;
                cCovs++;continue;
            }
                        if(cCovs==7 && itsCovs)
            {
                covTemp.pos_mvc[0][2]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }
                        if(cCovs==8 && itsCovs)
            {
                covTemp.pos_mvc[1][0]=stod(value);
           ////     cout<<value<<endl;
                cCovs++;continue;
            }
                                    if(cCovs==9 && itsCovs)
            {
                covTemp.pos_mvc[1][1]=stod(value);
           ////     cout<<value<<endl;
                cCovs++;continue;
            }
                                    if(cCovs==10 && itsCovs)
            {
                covTemp.pos_mvc[1][2]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }
                                if(cCovs==11 && itsCovs)
            {
                covTemp.pos_mvc[2][0]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }

                            if(cCovs==12 && itsCovs)
            {
                covTemp.pos_mvc[2][1]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }
                                        if(cCovs==13 && itsCovs)
            {
                covTemp.pos_mvc[2][2]=stod(value);
             ////   cout<<value<<endl;
                cCovs++;continue;
            }

                        if(cCovs==14 && itsCovs)
            {
                covTemp.vel_mvc[0][0]=stod(value);
             ////   cout<<value<<endl;
                cCovs++;continue;
            }
                        if(cCovs==15 && itsCovs)
            {
                covTemp.vel_mvc[0][1]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }
                        if(cCovs==16 && itsCovs)
            {
                covTemp.vel_mvc[0][2]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }
                        if(cCovs==17 && itsCovs)
            {
                covTemp.vel_mvc[1][0]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }
                                    if(cCovs==18 && itsCovs)
            {
                covTemp.vel_mvc[1][1]=stod(value);
             ////   cout<<value<<endl;
                cCovs++;continue;
            }
                                    if(cCovs==19 && itsCovs)
            {
                covTemp.vel_mvc[1][2]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }
                                if(cCovs==20 && itsCovs)
            {
                covTemp.vel_mvc[2][0]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }

                            if(cCovs==21 && itsCovs)
            {
                covTemp.vel_mvc[2][1]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }
                                        if(cCovs==22 && itsCovs)
            {
                covTemp.vel_mvc[2][2]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }

                        if(cCovs==23 && itsCovs)
            {
                covTemp.att_mvc[0][0]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }
                        if(cCovs==24 && itsCovs)
            {
                covTemp.att_mvc[0][1]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }
                        if(cCovs==25 && itsCovs)
            {
                covTemp.att_mvc[0][2]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }
                        if(cCovs==26 && itsCovs)
            {
                covTemp.att_mvc[1][0]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }
                                    if(cCovs==27 && itsCovs)
            {
                covTemp.att_mvc[1][1]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }
                                    if(cCovs==28 && itsCovs)
            {
                covTemp.att_mvc[1][2]=stod(value);
             ////   cout<<value<<endl;
                cCovs++;continue;
            }
                                if(cCovs==29 && itsCovs)
            {
                covTemp.att_mvc[2][0]=stod(value);
             ////   cout<<value<<endl;
                cCovs++;continue;
            }

                            if(cCovs==30 && itsCovs)
            {
                covTemp.att_mvc[2][1]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }
                                        if(cCovs==31 && itsCovs)
            {
                covTemp.att_mvc[2][2]=stod(value);
            ////    cout<<value<<endl;
                cCovs++;continue;
            }
            if(cCovs==32 && itsCovs)
            {
            //cov.push_back(covTemp);
            ////cout<<endl<<endl<<"FIM DE COVARIANCIAS"<<endl;
            cCovs=0;itsCovs=false;
            }

        }

}

void survData::photoTimestampReader()
{
// if the time in UTC format are without Leap Seconds, insert "0" as the total of leap seconds
//variables for time conversion
    unsigned short  tempWeek;
    double          tempTow;
    bool            convOK;

    //the log with the final timestamps
    ofstream finalTimes("momento_da_tomada.txt");

    //first: take away the separators
    //first: ":"
    ifstream rawlog(log2_);
    string line1;
    ofstream part1("p1.txt");
        while(getline(rawlog,line1,':'))
        {
        part1<<line1<<endl;
        }
        rawlog.close();part1.close();

        //first: ";"
        ifstream log2("p1.txt");
        string line2;
        ofstream part2("p2.txt");
        while(getline(log2,line2,';'))
        {
        //if(line2 == "Start" || line2 == "End"){part2<<line2;}
        //else {part2<<line2<<endl;}
        //if(line2 == "") cout<<"existe"<<endl;
        part2<<line2<<endl;
        }
        log2.close();part2.close();

    //now, we could read the one-value-per-line file
    ifstream log3("p2.txt");
    string line3;
    moment tempLS,tempLE,tempRS,tempRE;
    stereopairTimestamp tempStamp;
    double tHour,tMin,tSec;

    //control variables
    int count=1;
    bool OneLS=false;
    bool OneLE=false;
    bool OneRS=false;
    bool OneRE=false;
    //bool isLeft=false;
    //bool isRight=false;
    bool isLeft,isRight;
    int gravou1=0;

while(getline(log3,line3))
{
//if (line3 == "")cout<<line3<<endl;

        if (count==1)
        {
                if(line3 != "")
                {tHour=stod(line3);}
        count++;continue;
        }
        if (count==2)
        {
                if(line3 != "")
                {tMin=stod(line3);}
        count++;continue;
        }
        if (count==3)
        {
                if(line3 != "")
                {tSec=stod(line3);}
        count++;continue;
        }
        if (count==4)
        {
            if(line3=="1")
            {isLeft=true;isRight=false;}
            if(line3=="2")
            {isLeft=false;isRight=true;}
        count++;continue;
        }
        if (count==5)
        {
            if(line3=="Start")
            {
                if(isLeft)
                {
                    tempLS.hour=tHour;
                    tempLS.minute=tMin;
                    tempLS.second=tSec;
//                    cout<<"passo OK"<<endl;
                    if(!OneLS){OneLS=true;}
                }
                if(isRight)
                {
                    tempRS.hour=tHour;
                    tempRS.minute=tMin;
                    tempRS.second=tSec;
                    if(!OneRS){OneRS=true;}
                }
            }
            if(line3=="End")
            {
                if(isLeft && OneLS)
                {
                    tempLE.hour=tHour;
                    tempLE.minute=tMin;
                    tempLE.second=tSec;
                    if(!OneLE){OneLE=true;}
                }
                if(isRight && OneLS)
                {
                    tempRE.hour=tHour;
                    tempRE.minute=tMin;
                    tempRE.second=tSec;
                    if(!OneRE){OneRE=true;}
                }
            }
        count++;continue;
        }
        if (count==6)
        {
                if(OneLE && OneLS && OneRE && OneRS)
                {
                //record the timestamps in UTC
                tempStamp.leftStart     = tempLS;
                tempStamp.leftEnd       = tempLE;
                tempStamp.RightStart    = tempRS;
                tempStamp.RightEnd      = tempRE;
                tempStamp.tsyear        = year;
                tempStamp.tsmonth       = month;
                tempStamp.tsday         = day;
                tempStamp.leapSeconds=LeapSeconds;
                //record the timestamps in GPS time:
                convOK =  TIMECONV_GetGPSTimeFromRinexTime(year,month,day,
                tempLS.hour+float(toUTC),tempLS.minute,tempLS.second,&tempWeek,&tempTow);
                tempStamp._leftStart.week           =   tempWeek;
                tempStamp._leftStart.seconds        =   tempTow+LeapSeconds;

                convOK =  TIMECONV_GetGPSTimeFromRinexTime(year,month,day,
                tempLE.hour+float(toUTC),tempLE.minute,tempLE.second,&tempWeek,&tempTow);
                tempStamp._leftEnd.week              =  tempWeek;
                tempStamp._leftEnd.seconds           =   tempTow+LeapSeconds;

                convOK =  TIMECONV_GetGPSTimeFromRinexTime(year,month,day,
                tempRS.hour+float(toUTC),tempRS.minute,tempRS.second,&tempWeek,&tempTow);
                tempStamp._RightStart.week           = tempWeek;
                tempStamp._RightStart.seconds        = tempTow+LeapSeconds;

                convOK =  TIMECONV_GetGPSTimeFromRinexTime(year,month,day,
                tempRE.hour+float(toUTC),tempRE.minute,tempRE.second,&tempWeek,&tempTow);
                tempStamp._RightEnd.week            = tempWeek;
                tempStamp._RightEnd.seconds         = tempTow+LeapSeconds;

                //the final timestamp:
                tempStamp.finalTime.week            =tempWeek;
                tempStamp.finalTime.seconds         =( ((tempStamp._leftStart.seconds+tempStamp._leftEnd.seconds)/2) +  ((tempStamp._RightStart.seconds+tempStamp._RightEnd.seconds)/2) ) / 2;



                timestamps.push_back(tempStamp);
                OneLS=false;OneLE=false;OneRS=false;OneRE=false;
                gravou1++;
//                cout<<gravou1<<endl;
                }
        count=1;continue;
        }

} //end of the main process of the function

finalTimes << "semana GPS: " << tempWeek<<endl<<endl;


//cout<<timestamps.size()<<endl<<endl;
for (unsigned int i=0;i < timestamps.size();i++)
{
//cout<<endl<<i<<endl<<endl;
//cout<<timestamps.at(i).leftStart.hour<<endl;
//cout<<timestamps.at(i).leftStart.minute<<endl;
//cout<<timestamps.at(i).leftStart.second<<endl;
//
//cout<<timestamps.at(i).RightStart.hour<<endl;
//cout<<timestamps.at(i).RightStart.minute<<endl;
//cout<<timestamps.at(i).RightStart.second<<endl;
//
//cout<<timestamps.at(i).leftEnd.hour<<endl;
//cout<<timestamps.at(i).leftEnd.minute<<endl;
//cout<<timestamps.at(i).leftEnd.second<<endl;
//
//cout<<timestamps.at(i).RightEnd.hour<<endl;
//cout<<timestamps.at(i).RightEnd.minute<<endl;
//cout<<timestamps.at(i).RightEnd.second<<endl;
//
//cout<<timestamps.at(i).finalTime.week<<endl;
//cout<<std::setprecision(20)<<timestamps.at(i).finalTime.seconds<<endl<<endl;
//
//cout<<std::setprecision(20)<<timestamps.at(i)._leftStart.seconds<<endl;
//
//cout<<std::setprecision(20)<<timestamps.at(i)._leftEnd.seconds<<endl;
//
//cout<<std::setprecision(20)<<timestamps.at(i)._RightStart.seconds<<endl;
//
//cout<<std::setprecision(20)<<timestamps.at(i)._RightEnd.seconds<<endl<<endl;

//cout<<day<<"    "<<month<<"     "<<year<<endl;

finalTimes<<std::setprecision(13)<<timestamps.at(i).finalTime.seconds<<endl;
}

finalTimes.close();
}


// function for the interpolation of data
void survData::dataInterpolation()
{
double delta;
int stereopair = 0;
intervals tempIntervals;
//the differences in time:
ofstream deltas("diftempo.txt");
//the final positions
ofstream finalObservations("observacoes.txt");
//arbitrary values
double minNeg = -1000; int posMinNeg;
double minPos =  1000; int posMinPos;
//temporary variables
obsWcov tempObs;
obs tempObsBef,tempObsAf;
double tBef,tAft,time1;


        for (unsigned int i=0;i < timestamps.size();i++)
        {
        stereopair++;
        deltas<<"Estereopar  "<<stereopair<<endl<<endl;
            for(unsigned  int j=0;j < pvaWcov.size();j++)
            {
            delta =pvaWcov.at(j).observation.obsTime.seconds - timestamps.at(i).finalTime.seconds;

            if (delta <= 0 && delta > minNeg)
            {
            minNeg = delta;posMinNeg=j;
            }

            if (delta >= 0 && delta < minPos)
            {
            minPos = delta;posMinPos=j;
            }


            deltas<<std::setprecision(15)<<delta<<"  "<<minNeg<<"  "<<minPos<<endl;
            }
            deltas<<endl;

            deltas<<"menor negativo:  "<<minNeg<<"  menor positivo:  "<<minPos<<" primeira: "<<posMinNeg<<" segunda "<<posMinPos<<endl;
            deltas<<endl<<endl;
            //storing that shortest intervals and these positions
            tempIntervals.shortestNeg = minNeg;tempIntervals.shortestPos = minPos;
            tempIntervals.shortestNegPos = posMinNeg; tempIntervals.shortestPosPos = posMinPos;
            tempIntervals.theTime = timestamps.at(i).finalTime.seconds;
            finalIntervals.push_back(tempIntervals);
            //reassign the temporary values:
            minNeg = -1000;minPos=1000;
        }
        deltas.close();
        //now we have the intervals of time and the positions to do the interpolation

        for(unsigned int k=0;k < finalIntervals.size();k++)
        {
        tBef = pvaWcov.at(finalIntervals.at(k).shortestNegPos).observation.obsTime.seconds;
        tAft = pvaWcov.at(finalIntervals.at(k).shortestPosPos).observation.obsTime.seconds;
        time1 = finalIntervals.at(k).theTime;

        tempObsBef  =   pvaWcov.at(finalIntervals.at(k).shortestNegPos).observation;
        tempObsAf   =   pvaWcov.at(finalIntervals.at(k).shortestPosPos).observation;

        tempObs.observation.lat     = linTerp(time1,tBef,tAft,tempObsBef.lat,tempObsAf.lat);
        tempObs.observation.lgt     = linTerp(time1,tBef,tAft,tempObsBef.lgt,tempObsAf.lgt);
        tempObs.observation.h       = linTerp(time1,tBef,tAft,tempObsBef.h,tempObsAf.h);
        tempObs.observation.roll   = linTerp(time1,tBef,tAft,tempObsBef.roll,tempObsAf.roll);
        tempObs.observation.pitch     = linTerp(time1,tBef,tAft,tempObsBef.pitch,tempObsAf.pitch);
        tempObs.observation.azimuth   = linTerp(time1,tBef,tAft,tempObsBef.azimuth,tempObsAf.azimuth);
        tempObs.observation.velX    = linTerp(time1,tBef,tAft,tempObsBef.velX,tempObsAf.velX);
        tempObs.observation.velY    = linTerp(time1,tBef,tAft,tempObsBef.velY,tempObsAf.velY);
        tempObs.observation.velZ    = linTerp(time1,tBef,tAft,tempObsBef.velZ,tempObsAf.velZ);


        finalObservations.precision(40);
        finalObservations<<tempObs.observation.lat<<",";
        finalObservations<<tempObs.observation.lgt<<",";
        finalObservations<<tempObs.observation.h<<",";
        finalObservations<<tempObs.observation.roll<<",";
        finalObservations<<tempObs.observation.pitch<<",";
        finalObservations<<tempObs.observation.azimuth<<",";
        finalObservations<<tempObs.observation.velX<<",";
        finalObservations<<tempObs.observation.velY<<",";
        finalObservations<<tempObs.observation.velZ<<endl;

        //for the covariances, the NearestNeighbour is good enough
        if (fabs(finalIntervals.at(k).shortestNeg) < fabs(finalIntervals.at(k).shortestPos)){
        tempObs.obsCovariances = pvaWcov.at(finalIntervals.at(k).shortestNegPos).obsCovariances;}
        else {tempObs.obsCovariances = pvaWcov.at(finalIntervals.at(k).shortestPosPos).obsCovariances;}

        finalObsWcovs.push_back(tempObs);
        }
//cout<<finalObsWcovs.size()<<endl;



finalObservations.close();
}

