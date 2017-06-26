//#include <iostream>
//#include <ostream>
//#include <fstream>
//#include <armadillo>
//#include "master.hpp"
//
//Rzxy(,)
//F_A: functions for the Agisoft conventions
//double yaw_rotZXY(mat Rzxy,bool giveRadians=false)
//{
//
//    double c = (180/datum::pi);
//
//    if(giveRadians)
//    {
//        c = 1.0f;
//    }
//    F_A
//    return atan2(Rzxy(0,1),Rzxy(1,1)) * c;
//}
//
//double pitch_rotZXY(mat Rzxy,bool giveRadians=false)
//{
//    double c = (180/datum::pi);
//
//    if(giveRadians)
//    {
//        c = 1.0f;
//    }
//    F_A
//    return -asin(Rzxy(2,1)) * c;
//}
//
//double roll_rotZXY(mat Rzxy,bool giveRadians=false)
//{
//    double c = (180/datum::pi);
//
//    if(giveRadians)
//    {
//        c = 1.0f;
//    }
//    F_A
//    return atan2(Rzxy(2,0),Rzxy(2,2)) * c;
//}
//
//mat one_axis_rot(double a,int axis,bool giveInverseT=false,bool isRadians=false)
//{
//
//    these will give you the rotation in a CLOCKWISE way
//
//    mat res = eye(3,3);
//    double c = (datum::pi/180);
//
//    if (isRadians)
//    {
//        c = 1.0f;
//    }
//
//    double ca = cos(a*c);
//    double sa = sin(a*c);
//
//    if(axis == 1)
//    {
//        res(1,1) =  ca;
//        res(1,2) =  sa;
//        res(2,1) = -sa;
//        res(2,2) =  ca;
//    }
//    if(axis == 2)
//    {
//        res(0,0) =  ca;
//        res(0,2) = -sa;
//        res(2,0) =  sa;
//        res(2,2) =  ca;
//    }
//    if (axis == 3)
//    {
//        res(0,0) =  ca;
//        res(0,1) =  sa;
//        res(1,0) = -sa;
//        res(1,1) =  ca;
//    }
//
//    if (giveInverseT)
//    {res = res.t();}
//
//    return res;
//}
//
//mat novatel_DCM(double roll,double pitch,double azimuth,bool isRadians=false,bool giveInverseT=false)
//{
//    matrix that gives the transformation
//    FROM   BODY FRAME
//    TO     LOCAL LEVEL FRAME
//
//    accordingly to the manual, is a combination of the following order:
//    Rz * Rx * Ry  or R3 * R1 * R2 , all of them TRANSPOSED
//
//    NOVATEL definitions:
//    roll   is around the Y axis
//    pitch  is around the X axis
//    yaw    is around the Z axis (and also Azimuth)
//
//    there is the fourth argument that can be used to give the inverse transformation
//
//    mat res = zeros(3,3);
//
//    double c = (datum::pi/180);
//
//    if (isRadians)
//    {
//        c = 1.0f;
//    }
//
//    conversion between azimuth and yaw
//    double yaw = -azimuth;
//
//    double cr = cos(roll * c);
//    double sr = sin(roll * c);
//
//    double cp = cos(pitch * c);
//    double sp = sin(pitch * c);
//
//    double cy = cos(yaw * c);
//    double sy = sin(yaw * c);
//
//    linewise splitted
//    res(0,0)= cy*cr-sy*sp*sr;
//    res(0,1)= -sy*cp;
//    res(0,2)= cy*sr + sy*sp*cr;
//
//    res(1,0)=  sy*cr+cy*sp*sr;
//    res(1,1)=  cy*cp;
//    res(1,2)=  sy*sr-cy*sp*cr;
//
//    res(2,0)=-cp*sr;
//    res(2,1)= sp;
//    res(2,2)= cp*cr;
//
//    if (giveInverseT)
//    {res = res.t();}
//
//    return res;
//}
//
//
//struct outputPoseToPhotoscan
//{
//    struct for individual data
//    string img_name;
//    double lat,lgt,h,yaw,pitch,roll;
//    double lat0,lgt0,h0; //coordinates of the origin
//
//    vec3 v0ECEF;
//
//    mat attMat;
//    outputPoseToPhotoscan(vec3 vecBF,mat bsight,unsigned int ind);
//};
//
//outputPoseToPhotoscan::outputPoseToPhotoscan(vec3 vBF,mat bsight,unsigned int ind)
//{
//    constants for multiplication and transformation
//    double to_deg = 180/datum::pi;
//    double to_rad = datum::pi/180;
//
/// First: position of the camera CP
//    filling the origin
//    lat0 = job.finalObsWcovs.at(ind).observation.lat;
//    lgt0 = job.finalObsWcovs.at(ind).observation.lgt;
//    h0   = job.finalObsWcovs.at(ind).observation.h;
//
//    converting the origin to XYZ
//    GEODESY_ConvertGeodeticCurvilinearToEarthFixedCartesianCoordinates(
//    GEODESY_REFERENCE_ELLIPSE_WGS84,
//    lat0*to_rad,lgt0*to_rad,h0,&v0ECEF(0),&v0ECEF(1),&v0ECEF(2));
//
//    compute the rotation matrix
//    mat R = novatel_DCM(
//    job.finalObsWcovs.at(ind).observation.roll,
//    job.finalObsWcovs.at(ind).observation.pitch,
//    job.finalObsWcovs.at(ind).observation.azimuth);
//
//    transform the vector to the LLF
//    vec3 vLLF = R * vBF;
//    transform to the ECEF system
//    vec3 vecECEF = ENU_to_ECEF(vLLF,v0ECEF,lat0,lgt0,h0);
//
//    transform to lat, long, h
//    GEODESY_ConvertEarthFixedCartesianToGeodeticCurvilinearCoordinates(
//    GEODESY_REFERENCE_ELLIPSE_WGS84,vecECEF(0),vecECEF(1),vecECEF(2),&lat,&lgt,&h);
//
//    converting to degrees:
//    lat *= to_deg;lgt *= to_deg;
//
//    / second: orientation
//
//         R is from IMU BF to IMU CN LLF
//
//        matrix from ECEF to IMU CN LLF
//        mat Rel1 = R_ecef_enu(lat0,lgt0);
//
//        matrix from ECEF to camera CP LLF
//        mat Rel2 = R_ecef_enu(lat,lgt);
//
//        matrix from IMU CN LLF to camera CP LLF
//        mat Rl1l2 = Rel1 * Rel2.t();
//
//        matrix from the camera CP LLF to IMU BF
//        mat Rl2bf1 = R.t() * Rl1l2.t();
//
//        finally, the camera BF to camera CP LLF
//        the bsight needs to be Rbf2bf1, aka from Camera BF to IMU
//        attMat = Rl2bf1.t()*bsight.t();
//
//        now the yaw pitch roll to photoscan
//        yaw     = yaw_rotZXY(attMat);
//        pitch   = pitch_rotZXY(attMat);
//        roll    = roll_rotZXY(attMat);
//    }
//
//struct outputterToPhotoscan
//{
// vector<outputPoseToPhotoscan> individuals;
//
//
//
// outputterToPhotoscan(void);
//};
//
//outputterToPhotoscan::outputterToPhotoscan(void)
//{
// ofstream outL("to_photoscan_left.txt");
// ofstream outR("to_photoscan_right.txt");
//
// string separator = "   ";
//
//    for (unsigned int i=0;i<job.finalObsWcovs.size();i++)
//    {
//        outputPoseToPhotoscan  Left(rCalib.LcamLA,rCalib.Rimu_LC.t(),i);
//        outputPoseToPhotoscan Right(rCalib.RcamLA,rCalib.Rimu_RC.t(),i);
//
//        outL<<Left.img_name<<separator<<Left.lat<<separator<<Left.lgt<<separator<<Left.h;
//        outL<<separator<<Left.yaw<<separator<<Left.pitch<<separator<<Left.roll<<endl;
//
//        outR<<Right.img_name<<separator<<Right.lat<<separator<<Right.lgt<<separator<<Left.h;
//        outR<<separator<<Right.yaw<<separator<<Right.pitch<<separator<<Right.roll<<endl;
//    }
//}
//
