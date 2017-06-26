#include <iostream>
#include <ostream>
#include <fstream>
//#include "logreader.hpp"
#include "logreader2.hpp"
#include "or_relativ.hpp"
#include "ecef_enu.hpp"
#include "geodesy.c"
#include "geodesy.h"
#include "photoreader.hpp"
#include "timeChecker.hpp"

//allways remember: the code follows the c++11 standard

using namespace std;

///------------***------------///

//doing the calibration for relative orientation
//SMMTleverARM rCalib("or/data.tcl");

//    SMMTleverARM rCalib("/home/kauevestena/Dropbox/IC/cut_010916/data2.tcl");
//
//    //reading the data from navigation
//    //survData job("nav/volta2.txt","nav/simulados.txt",2015,03,20,16);
//    survData job(
//    "/home/kauevestena/Dropbox/IC/cut_010916/nav.ASC",
//    "/home/kauevestena/Dropbox/IC/cut_010916/01_09_2016_SHOT.csv",2016,9,1,0);
//
//
//    photoReader jobPhotos("/home/kauevestena/Dropbox/IC/cut_010916","*.JPG",
//    "esquerda","direita","/home/kauevestena/Dropbox/IC/cut_010916/out",
//    "/home/kauevestena/Dropbox/IC/cut_010916/intrinsics.yml");
///------------***------------///

///MISC functions



string subsBefLChar(string input,string character,bool includeChar=false)
{
//return the string before the last occur. of the specified character
    int pos = input.rfind(character);

    if(pos != -1)
    {
        if (includeChar)
        {
            return input.substr(0,pos)+character;
        }
        else
        {
            return input.substr(0,pos);
        }

    }

    else
    {
        return "";
    }

}

string subsAfLChar(string input,string character)
{
//return the string after the last occur. of the specified character
    int pos = input.rfind(character);
    if(pos != -1)
    {
        return input.substr(pos+character.size(),string::npos);
    }
    else
    {
        return "";
    }
}
///end of MISC functions

//Rzxy(,)
//F_A: functions for the Agisoft conventions
double yaw_rotZXY(mat Rzxy,bool giveRadians=false)
{
    double c = (180/datum::pi);

    if(giveRadians)
    {
        c = 1.0f;
    }
    //F_A
//    return atan2(Rzxy(0,1),Rzxy(1,1)) * c;

    //gimbal lock caring

////    if (Rzxy(2,1) > 0.999 || Rzxy(2,1) < -0.999)
////    {
////        return atan2(Rzxy(1,0),Rzxy(0,0)) * c;
////    }

    double res = atan2(Rzxy(0,1),Rzxy(1,1)) * c;

    if (res < 0)
    {
        res+=360;
    }

    return res;
}

double pitch_rotZXY(mat Rzxy,bool giveRadians=false)
{
    double c = (180/datum::pi);

    if(giveRadians)
    {
        c = 1.0f;
    }
    //F_A
    //return -asin(Rzxy(2,1)) * c;

    //gimbal lock caring

////    if (Rzxy(2,1) > 0.999)
////    {
////        return datum::pi/2 * c;
////    }
////
////    if (Rzxy(2,1) < -0.999)
////    {
////        return -datum::pi/2 * c;
////    }

    return asin(Rzxy(2,1)) * c;
}

double roll_rotZXY(mat Rzxy,bool giveRadians=false)
{
    double c = (180/datum::pi);

    if(giveRadians)
    {
        c = 1.0f;
    }
    //F_A
//    return -atan2(Rzxy(2,0),Rzxy(2,2)) * c;

    //gimbal lock caring


////    if (Rzxy(2,1) > 0.999 || Rzxy(2,1) < -0.999)
////    {
////        return 0;
////    }


    return atan2(-Rzxy(2,0),Rzxy(2,2)) * c;
}

//for PIX4D convention
//omega = atan2(-M(2, 3),M(3, 3))
//phi   = atan2(M(1, 3),sqrt(M(2, 3)² + M(3, 3)²))
//phi   = asin(M(1, 3))
//kappa = atan2(-M(1, 2),M(1, 1))

double omega_rotXYZ(mat Rxyz,bool giveRadians=false)
{
    double c = (180/datum::pi);

    if(giveRadians)
    {
        c = 1.0f;
    }

    return atan2(-Rxyz(1,2),Rxyz(2,2)) * c;
}

double phi_rotXYZ(mat Rxyz,bool giveRadians=false)
{
    double c = (180/datum::pi);

    if(giveRadians)
    {
        c = 1.0f;
    }

    return asin(Rxyz(0,2)) * c;
}

double kappa_rotXYZ(mat Rxyz,bool giveRadians=false)
{
    double c = (180/datum::pi);

    if(giveRadians)
    {
        c = 1.0f;
    }

    return atan2(-Rxyz(0,1),Rxyz(0,0)) * c;
}

mat one_axis_rot(double a,int axis,bool giveInverseT=false,bool isRadians=false)
{

    //these will give you the rotation in a CLOCKWISE way

    mat res = eye(3,3);
    double c = (datum::pi/180);

    if (isRadians)
    {
        c = 1.0f;
    }

    double ca = cos(a*c);
    double sa = sin(a*c);

    if(axis == 1)
    {
        res(1,1) =  ca;
        res(1,2) =  sa;
        res(2,1) = -sa;
        res(2,2) =  ca;
    }
    if(axis == 2)
    {
        res(0,0) =  ca;
        res(0,2) = -sa;
        res(2,0) =  sa;
        res(2,2) =  ca;
    }
    if (axis == 3)
    {
        res(0,0) =  ca;
        res(0,1) =  sa;
        res(1,0) = -sa;
        res(1,1) =  ca;
    }

    if (giveInverseT)
    {
        res = res.t();
    }

    return res;
}

mat novatel_DCM(double roll,double pitch,double azimuth,bool isRadians=false,bool giveInverseT=false)
{
    //matrix that gives the transformation
    //FROM   BODY FRAME
    //TO     LOCAL LEVEL FRAME

    //accordingly to the manual, is a combination of the following order:
    //Rz * Rx * Ry  or R3 * R1 * R2 , all of them TRANSPOSED

    //NOVATEL definitions:
    //roll   is around the Y axis
    //pitch  is around the X axis
    //yaw    is around the Z axis (and also Azimuth)

    //there is the fourth argument that can be used to give the inverse transformation

    mat res = zeros(3,3);

    double c = (datum::pi/180);

    if (isRadians)
    {
        c = 1.0f;
    }

    //conversion between azimuth and yaw
//    double yaw = -azimuth;
    double yaw = 360-azimuth;

    double cr = cos(roll * c);
    double sr = sin(roll * c);

    double cp = cos(pitch * c);
    double sp = sin(pitch * c);

    double cy = cos(yaw * c);
    double sy = sin(yaw * c);

    //linewise splitted
    res(0,0)= cy*cr-sy*sp*sr;
    res(0,1)= -sy*cp;
    res(0,2)= cy*sr + sy*sp*cr;

    res(1,0)=  sy*cr+cy*sp*sr;
    res(1,1)=  cy*cp;
    res(1,2)=  sy*sr-cy*sp*cr;

    res(2,0)=-cp*sr;
    res(2,1)= sp;
    res(2,2)= cp*cr;

    if (giveInverseT)
    {
        res = res.t();
    }

    return res;
}


struct outputPoseToPhotoscan
{
    //struct for individual data
    //string img_name;
    double lat,lgt,h,yaw,pitch,roll,omega_p4d,phi_p4d,kappa_p4d;
    double lat0,lgt0,h0; //coordinates of the origin

    vec3 v0ECEF;

    mat attMat;
    outputPoseToPhotoscan(vec3 vecBF,mat bsight,unsigned int ind,survData2 job);
};

outputPoseToPhotoscan::outputPoseToPhotoscan(vec3 vBF,mat bsight,unsigned int ind,survData2 job)
{
    //constants for multiplication and transformation
    double to_deg = 180/datum::pi;
    double to_rad = datum::pi/180;

/// First: position of the camera CP
    //filling the origin
    lat0 = job.pvaWcov.at(ind).observation.lat;
    lgt0 = job.pvaWcov.at(ind).observation.lgt;
    h0   = job.pvaWcov.at(ind).observation.h;

    //converting the origin to XYZ
    GEODESY_ConvertGeodeticCurvilinearToEarthFixedCartesianCoordinates(
        GEODESY_REFERENCE_ELLIPSE_WGS84,
        lat0*to_rad,lgt0*to_rad,h0,&v0ECEF(0),&v0ECEF(1),&v0ECEF(2));

    //compute the rotation matrix
    mat R = novatel_DCM(
                job.pvaWcov.at(ind).observation.roll,
                job.pvaWcov.at(ind).observation.pitch,
                job.pvaWcov.at(ind).observation.azimuth);

    //transform the vector to the LLF
    vec3 vLLF = R * vBF;

    //cout<<endl<<vLLF<<endl;

    //transform to the ECEF system
    vec3 vecECEF = ENU_to_ECEF(vLLF,v0ECEF,lat0,lgt0);

    //transform to lat, long, h
    GEODESY_ConvertEarthFixedCartesianToGeodeticCurvilinearCoordinates(
        GEODESY_REFERENCE_ELLIPSE_WGS84,vecECEF(0),vecECEF(1),vecECEF(2),&lat,&lgt,&h);

    //converting to degrees:
    lat *= to_deg;
    lgt *= to_deg;

    /// second: orientation

    // R is from IMU BF to IMU CN LLF

    //matrix from ECEF to IMU CN LLF
    mat Rel1 = R_ecef_enu(lat0,lgt0);

    //matrix from ECEF to camera CP LLF
    mat Rel2 = R_ecef_enu(lat,lgt);

    //matrix from IMU CN LLF to camera CP LLF
    mat Rl1l2 = Rel1 * Rel2.t();

    //matrix from the camera CP LLF to IMU BF
    mat Rl2bf1 = R.t() * Rl1l2.t();

    //finally, the camera BF to camera CP LLF
    //the bsight needs to be Rbf2bf1, aka from Camera BF to IMU
    attMat = Rl2bf1.t()*bsight.t();

    //now the yaw pitch roll to photoscan
    yaw     = yaw_rotZXY(attMat);
    pitch   = pitch_rotZXY(attMat);
    roll    = roll_rotZXY(attMat);

    //omega phi kappa for pix4d
    omega_p4d = omega_rotXYZ(attMat);
    phi_p4d = phi_rotXYZ(attMat);
    kappa_p4d = kappa_rotXYZ(attMat);
}

struct outputterToPhotoscan
{
    vector<outputPoseToPhotoscan> individuals;

//     SMMTleverARM rCalib;
//
//    survData job;
//
//    photoReader jobPhotos;


    outputterToPhotoscan(SMMTleverARM rCalib,survData2 job,photoReader jobPhotos,bool sameConvention,bool oldConvention);
};

outputterToPhotoscan::outputterToPhotoscan(SMMTleverARM rCalib,survData2 job,photoReader jobPhotos,bool sameConvention = false,bool oldConvention = false)
{
    double to_deg = 180/datum::pi;
    double to_rad = datum::pi/180;

    string separator = ",";

    //cout<<"reached 02"<<endl;

    string outpathL = subsBefLChar(jobPhotos.outNamesL.at(0),"/",true)+"to_photoscan_left.txt";
    string outpathR = subsBefLChar(jobPhotos.outNamesR.at(0),"/",true)+"to_photoscan_right.txt";

    string outpathL_enu =  subsBefLChar(jobPhotos.outNamesL.at(0),"/",true)+"enu_photoscan_L.txt";
    string outpathR_enu =  subsBefLChar(jobPhotos.outNamesR.at(0),"/",true)+"enu_photoscan_R.txt";

    string outpathL_pix4d =  subsBefLChar(jobPhotos.outNamesL.at(0),"/",true)+"to_pix4d_L.txt";
    string outpathR_pix4d =  subsBefLChar(jobPhotos.outNamesR.at(0),"/",true)+"to_pix4d_R.txt";

    cout<<outpathL<<endl;
    cout<<outpathR<<endl<<endl;
    cout<<outpathL_enu<<endl;
    cout<<outpathR_enu<<endl;

    ofstream outL(outpathL);
    outL.precision(25);
    ofstream outR(outpathR);
    outR.precision(25);

    ofstream outL_enu(outpathL_enu);
    outL_enu.precision(25);
    ofstream outR_enu(outpathR_enu);
    outR_enu.precision(25);

    ofstream outL_pix4D(outpathL_pix4d);
    outL_pix4D.precision(25);
    ofstream outR_pix4D(outpathR_pix4d);
    outR_pix4D.precision(25);


// cout<<outpathL<<endl<<outpathR<<endl<<endl;

    unsigned t_photos    =   jobPhotos.outNamesL.size();
    unsigned t_obs       =    job.pvaWcov.size();
//    bool tm_obs = false; //tm: too much
//    bool tm_photos = false;

    cout<<endl<<"total de obs: "<<t_obs<<", total de fotos: "<<t_photos<<endl<<endl;

    if (t_obs>t_photos)
    {
        cout<<"Ha "<<t_obs-t_photos<<" mais observacoes que fotos, sera cortado"<<endl;
//        tm_obs = true;
    }
    else if (t_obs<t_photos)
    {
        cout<<"Ha "<<t_photos-t_obs<<" menos observacoes que fotos, sera cortado"<<endl;
//        tm_photos = true;
    }

    unsigned int j;

    /// ***
    //establishing an origin for the ENU output
    vec3 P0_ECEF;

    double lat0 =   job.pvaWcov.at(0).observation.lat;
    double lgt0 =   job.pvaWcov.at(0).observation.lgt;
    double h0 =     job.pvaWcov.at(0).observation.h;

        GEODESY_ConvertGeodeticCurvilinearToEarthFixedCartesianCoordinates(
        GEODESY_REFERENCE_ELLIPSE_WGS84,
        lat0*to_rad,lgt0*to_rad,h0,&P0_ECEF(0),&P0_ECEF(1),&P0_ECEF(2));
    /// ***



    for (unsigned int i=0; i<job.pvaWcov.size(); i++)
    {


        //boresighting matrices
        mat BsightL = rCalib.Rimu_LC.t();
        mat BsightR = rCalib.Rimu_RC.t();

//    cout<<BsightL<<endl<<endl;

        j = i;

        if (sameConvention)
        {
//    cout<<"reached"<<endl<<endl;
            BsightL = rCalib.Rimu_LC2.t();
            BsightR = rCalib.Rimu_RC2.t();
        }

        if (oldConvention && !sameConvention)
        {
            BsightL = rCalib.Rimu_LC3.t();
            BsightR = rCalib.Rimu_RC3.t();

            cout<<endl<<"utilizando a convenção fotogrametrica clássica (eixo optico para trás)"<<endl;
        }

        if (!oldConvention && !sameConvention)
        {
            cout<<endl<<"utilizando a convenção fotogrametrica moderna (eixo optico para frente)"<<endl;
        }

         if (!oldConvention && sameConvention)
        {
            cout<<endl<<"utilizando mesma convenção da IMU (eixo Y para frente)"<<endl;
        }

//    cout<<BsightL<<endl<<endl;

        outputPoseToPhotoscan  Left(rCalib.LcamLA,BsightL,i,job);
        outputPoseToPhotoscan Right(rCalib.RcamLA,BsightR,i,job);

    /// FOR ENU OUTPUT
    //ECEF vectors
    vec3 ecefL,ecefR;

        GEODESY_ConvertGeodeticCurvilinearToEarthFixedCartesianCoordinates(
        GEODESY_REFERENCE_ELLIPSE_WGS84,
        Left.lat*to_rad,Left.lgt*to_rad,Left.h,&ecefL(0),&ecefL(1),&ecefL(2));

        GEODESY_ConvertGeodeticCurvilinearToEarthFixedCartesianCoordinates(
        GEODESY_REFERENCE_ELLIPSE_WGS84,
        Right.lat*to_rad,Right.lgt*to_rad,Right.h,&ecefR(0),&ecefR(1),&ecefR(2));

    //ENU vectors
    vec3 enuL = ECEF_to_ENU(ecefL,P0_ECEF,lat0,lgt0);
    vec3 enuR = ECEF_to_ENU(ecefR,P0_ECEF,lat0,lgt0);

    /// ***

//        cout<<subsAfLChar(jobPhotos.outNamesL.at(i),"/")<<separator<<Left.lat<<separator<<Left.lgt<<separator<<Left.h;
//        cout<<separator<<Left.yaw<<separator<<Left.pitch<<separator<<Left.roll<<endl<<endl;
//
//        cout<<subsAfLChar(jobPhotos.outNamesR.at(i),"/")<<separator<<Right.lat<<separator<<Right.lgt<<separator<<Left.h;
//        cout<<separator<<Right.yaw<<separator<<Right.pitch<<separator<<Right.roll<<endl<<endl;

        vec3 posStd;
        posStd(0) = job.pvaWcov.at(i).obsCovariances.pos_mvc[0][0];
        posStd(1) = job.pvaWcov.at(i).obsCovariances.pos_mvc[1][1];
        posStd(2) = job.pvaWcov.at(i).obsCovariances.pos_mvc[2][2];

        //the fakes again
        bool fakeLeft   = false;
        bool fakeRight  = false;

        if (jobPhotos.outNamesL.at(j).find("fake") < jobPhotos.outNamesL.at(j).size() || jobPhotos.outNamesL.at(j).find("fake") < 0)
        {
            fakeLeft    = true;
        }

        if (jobPhotos.outNamesR.at(j).find("fake") < jobPhotos.outNamesR.at(j).size() || jobPhotos.outNamesR.at(j).find("fake") < 0)
        {
            fakeRight   = true;
        }


        if (!fakeLeft)
        {
            //LLh output
            outL<<subsAfLChar(jobPhotos.outNamesL.at(j),"/")<<separator<<Left.lgt<<separator<<Left.lat<<separator<<Left.h;
            outL<<separator<<Left.yaw<<separator<<Left.pitch<<separator<<Left.roll<<separator<<norm(arma::sqrt(posStd))<<endl;

            //ENU output
            outL_enu<<subsAfLChar(jobPhotos.outNamesL.at(j),"/")<<separator<<enuL(0)<<separator<<enuL(1)<<separator<<enuL(2);
            outL_enu<<separator<<Left.yaw<<separator<<Left.pitch<<separator<<Left.roll<<separator<<norm(arma::sqrt(posStd))<<endl;

            //p4d output
            outL_pix4D<<subsAfLChar(jobPhotos.outNamesL.at(j),"/")<<separator<<Left.lat<<separator<<Left.lgt<<separator<<Left.h;
            outL_pix4D<<separator<<Left.omega_p4d<<separator<<Left.phi_p4d<<separator<<Left.kappa_p4d<<separator<<norm(arma::sqrt(posStd))<<endl;
        }

        if (!fakeRight)
        {
            //LLh output
            outR<<subsAfLChar(jobPhotos.outNamesR.at(j),"/")<<separator<<Right.lgt<<separator<<Right.lat<<separator<<Left.h;
            outR<<separator<<Right.yaw<<separator<<Right.pitch<<separator<<Right.roll<<separator<<norm(arma::sqrt(posStd))<<endl;

            //ENU output
            outR_enu<<subsAfLChar(jobPhotos.outNamesR.at(j),"/")<<separator<<enuR(0)<<separator<<enuR(1)<<separator<<enuR(2);
            outR_enu<<separator<<Right.yaw<<separator<<Right.pitch<<separator<<Right.roll<<separator<<norm(arma::sqrt(posStd))<<endl;

            //p4d output
            outR_pix4D<<subsAfLChar(jobPhotos.outNamesR.at(j),"/")<<separator<<Right.lat<<separator<<Right.lgt<<separator<<Right.h;
            outR_pix4D<<separator<<Right.omega_p4d<<separator<<Right.phi_p4d<<separator<<Right.kappa_p4d<<separator<<norm(arma::sqrt(posStd))<<endl;
        }

        if (i+1 == t_photos || i+1 == t_obs)
        {
            break;
        }

    }

    outL.close();
    outR.close();

}

int main()
{
//    survData job_(
//        "/home/kauevestena/Dropbox/IC/cut_010916/nav.ASC",
//        "/home/kauevestena/Dropbox/IC/cut_010916/01_09_2016_SHOT.csv",2016,9,1,0);

survData2 job_("/home/kauevestena/Link to Google Drive/dados_SMMT/lev06122016/lev06122016 _mod.txt",true);

//    survData2 job_("/home/kauevestena/Link to Google Drive/dados_SMMT/lev06122016/lev06122016 _mod.txt",true);

//

    SMMTleverARM rCalib_("/home/kauevestena/Dropbox/IC/current_data/data2.tcl");
//
//    photoReader jobPhotos_("/home/kauevestena/Link to Google Drive/dados_SMMT/lev01112016/","*.JPG",
//                           "esquerda","esquerda","/home/kauevestena/data/undistorted",
//                           "/home/kauevestena/Dropbox/IC/current_data/intrinsics.yml");
//
//        photoReader jobPhotos_("/home/kauevestena/data/lev01112016/","*.JPG",
//                           "esquerda","direita","/home/kauevestena/data/undistorted",
//                           "/home/kauevestena/Dropbox/IC/current_data/intrinsics.yml");

//    photoReader jobPhotos_("/home/kauevestena/Link to Google Drive/dados_SMMT/lev06122016/","*.JPG",
//                           "esquerda","direita","/home/kauevestena/data/undistorted",
//                           "/home/kauevestena/Dropbox/IC/current_data/intrinsics.yml");

   photoReader jobPhotos_("/home/kauevestena/Link to Google Drive/dados_SMMT/lev06122016/","*.JPG",
                           "esquerda","direita","/home/kauevestena/data/teste",
                           "/home/kauevestena/Dropbox/IC/current_data/intrinsics.yml",false);


//
    if (job_.bugstatus == false)
    {
        outputterToPhotoscan(rCalib_,job_,jobPhotos_);
    }
    else
    {
        return -1;
    }


    return 0;
}
