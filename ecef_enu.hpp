//functions to convert from ECEF to ENU and vice-versa
#include <armadillo>

using namespace arma;

mat R_ecef_enu(double lat,double lgt,bool give_enu_ecef=false,bool isRadians=false)
{
//in the direct version, we are working with the matrix that
//rotates from ECEF to ENU. If you want from ENU to ECEF, use the third parameter.
//font: https://en.m.wikipedia.org/wiki/Geographic_coordinate_conversion#Molodensky_transformation

    mat res = eye(3,3);

    double c = datum::pi/180;

    if(isRadians)
    {
    c= 1.0f;
    }

    double sLat = sin(lat*c),cLat = cos(lat*c);
    double sLgt = sin(lgt*c),cLgt = cos(lgt*c);

    res(0,0)=-sLgt;        res(0,1)= cLgt;         res(0,2)=0;
    res(1,0)=-sLat*cLgt;   res(1,1)=-sLat*sLgt;    res(1,2)= cLat;
    res(2,0)= cLat*cLgt;   res(2,1)= cLat*sLgt;    res(2,2)= sLat;

    if(give_enu_ecef)
    {
    res = res.t();
    }

    return res;
}

vec3 ECEF_to_ENU(vec3 vECEF,vec3 v0ECEF,double lat0,double lgt0,bool givenRadians = false)
{
    //as output an ENU vector
    mat R = R_ecef_enu(lat0,lgt0,false,givenRadians);
    //cout<<"ECEFtoENU"<<endl<<R<<endl;
    return R*(vECEF-v0ECEF);
}
//

vec3 ENU_to_ECEF(vec3 vENU,vec3 v0ECEF,double lat0,double lgt0,bool givenRadians = false)
{
    //as output an ECEF vector
    mat R = R_ecef_enu(lat0,lgt0,true,givenRadians);
    //cout<<"ENUtoECEF"<<endl<<R*R.t()<<endl;
    return (R*vENU)+v0ECEF;
}

/// end of the file
