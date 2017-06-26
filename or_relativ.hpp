#include <iostream>
#include <armadillo>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <iomanip>


using namespace std;
using namespace arma;

struct GMS
{
    int G,M,S;

    GMS(int,int,int);
    double gdec();
    double rad();
};

struct topoPoint
{
    double Hz,Ze,Di;
    vec3 XYZ;
    bool rad;

    topoPoint(double,double,double,bool);
    void cart();
    void translate(vec3 tvec);
    void Rotate(mat R);

};


GMS::GMS (int g,int m,int s)
{
    G=g;
    M=m;
    S=s;
}

double GMS::gdec ()
{
    return double(G)+ (double(M)/60) + (double(S)/3600);
}

double GMS::rad()
{
    double radd = gdec()  * (datum::pi/180);
    return radd;
}

struct SMMTleverARM
{
//to handle with the base from the file
    vector<vector<string>> data;
    vector<vector<GMS>>  hor,zen;
    vector<vector<double>> dists;

//the arbitrary-frame surveyed points
    vector<vector<topoPoint>> points;

//the planes of the process (coefficients a,b,c,d)
    vec4 Hplan,Vplan;
    vec4 CEin,CEout,CDin,CDout;
//aux1 is parallel to Hplan

//matrices containing points
    mat ptsH,ptsV; //at the sides of the IMU
    mat ptsCEin,ptsCEout,ptsCDin,ptsCDout; //points on the parallel circular patterns, on the lens of the cameras
    mat ptsOnCEin,ptsOnCEout,ptsOnCDin,ptsOnCDout; //adjusted points

//rotation matrices
    mat Rot1 = eye<mat>(3,3);

//the principal vectors, to define the cannonical system
    vec3 u,v,w; //u for x,v for y,w for z

//the axis direction of each camera system, in the same convention of the IMU
    vec3 xE,yE,zE,xD,yD,zD;

//the axis of each camera system, with the z axis pointing backwards
    vec3 Xe,Ye,Ze,Xd,Yd,Zd;

//the difference between the antena top and the center of phase
    vec3 PhC = {0,0,-0.0071};

//the lever-arm vectors
    vec3 antLA,LcamLA,RcamLA;

//temporary, and auxiliar vectors
    vec3 tempE,tempD,aux1e,aux1d,aux_E0,aux_D0;

//the boresight matrices
    mat Rimu_LC = eye(3,3),Rimu_RC = eye(3,3),Rimu_LC2 = eye(3,3),Rimu_RC2 = eye(3,3),Rimu_LC3 = eye(3,3),Rimu_RC3 = eye(3,3);

//center of the IMU, the mark on the H plane, the mark on the v plane
    vec3 imuC,Hm,Vm;


    double planesAngle;

    SMMTleverARM (string fileName);


    vector<vector<string>> readLAFile(string filename);

//the member functions
    void horiz();
    void zenit();
    void distanc();
    void avgOBS();
    void pointsPrint(string filename);
    void translateAll(vec3 tvec);
    void rotAll (mat R);
    void report();

};

vec3 onPlanePoint(vec4 plCoef,double x,double y);

//all the function propotypes:
double pdpiH(double,double);
double pdpiV(double,double);
void TptMat (mat *M,vec Tvec);
void RptMat (mat *M,mat R);

vec4 plane3points(vec3 u,vec3 v,vec3 w);

mat topoPointsMat (vector<topoPoint> points,uword first,uword last);

vec4 parallelPlane(vec4 plane1,vec3 pointOnPlane);

vec3 ProjPtOrtPlane(vec3 pt2Proj,vec3 plNorm,vec3 ptOnPlane);

vec4 leasqPlane2(vector<topoPoint> pointList,mat *ptOnPlane,string repNam,uword ijpt=0)
{
    vec4 res,Xo;
    mat data,P,Pinf,A,B,M,mvcXa,mvcLa,mvcV,mvcW;
    vec Lb,W,X,Xa,K,V,La;
    double vp;

    if (pointList.size() >= 3)
    {
        Xo = plane3points(pointList.at(0).XYZ,pointList.at(1).XYZ,pointList.at(2).XYZ);
    }


    Pinf = eye(3,3)*10000;

    if (pointList.size() < 3)
    {
        res = {0,0,1,0};
    }
    else if (pointList.size() == 3)
    {
        res = Xo;
    }
    else
    {
        ofstream out(repNam);


        //vetor Lb
        data = ones(pointList.size(),3);

        for (uword i0 = 0; i0 < pointList.size(); i0++)
        {
            data.row(i0) = pointList.at(i0).XYZ.t();
        }

        Lb = vectorise(data.t());

        //matriz dos pesos
        P = eye(Lb.n_elem,Lb.n_elem);
        if (ijpt != 0)
        {
            //em caso de que haja um ponto a ser injuncionado
//        P.submat(ijpt*3-2,ijpt*3-2,ijpt*3,ijpt*3) = Pinf;
            P.submat(ijpt*3-3,ijpt*3-3,ijpt*3-1,ijpt*3-1) = Pinf;
        }
        //matriz jacobiana A (parametros)

        A = ones<mat>(pointList.size(),4);

        for (uword i = 0; i < pointList.size(); i++)
        {
            A(i,0)=data(i,0);
            A(i,1)=data(i,1);
            A(i,2)=data(i,2);
        }

        //matriz jacobiana B (observacoes)

        B = zeros<mat>(pointList.size(),Lb.n_elem);

        uword j = 0;
        for (uword i2 = 0; i2 < pointList.size(); i2++)
        {
            B(i2,0+j)=Xo(0);
            B(i2,1+j)=Xo(1);
            B(i2,2+j)=Xo(2);

            j += 3;
        }

        W = zeros<vec>(pointList.size());

        for (uword i3 = 0; i3 < pointList.size(); i3++)
        {
            W(i3) = dot(A.row(i3),Xo);
        }

        M = B * P.i() * B.t();

        X = - inv(A.t()*M.i()*A)*(A.t()*M.i()*W);

        Xa = Xo + X;

        K = -M.i()*(A*X+W);

        V = P.i() * B.t() * K;

        La = Lb+V;

        *ptOnPlane = trans(reshape(La,3,pointList.size()));

        vp = as_scalar( V.t() * P * V / (Lb.n_elem - Xo.n_elem) );

        mvcXa = vp * inv(A.t()*M.i()*A);

        mvcLa = vp*(inv(P)+inv(P)*B.t()*inv(M)*A*inv(A.t()*inv(M)*A)*A.t()*inv(M)*B*inv(P)-inv(P)*B.t()*inv(M)*B*inv(P));

        mvcV  = vp*P.i()-mvcLa;

        mvcW  = vp*M;

        out.precision(30);


        res = Xa;

//        Xo.raw_print(out,"Xo");
//        X.raw_print(out,"X");
//        Xa.raw_print(out,"Xa:");
//        res.raw_print(out," b ");

//        out <<"Xo+X: "<<endl<<std::setprecision(30)<< Xo+X<<endl;
        //normalizing the plane equations
        res.rows(0,2) = normalise(Xa.rows(0,2));
        res(3) = - dot(normalise(Xa.rows(0,2)),La.rows(0,2));

//            res.print(out," a ");

        data.print(out,"dados de entrada: ");
        out<<endl;

        res.print(out,"coeficientes do plano estimado: ");
        out<<endl;

        V.print(out,"resíduos: ");
        out<<endl;

        La.print(out,"La:");
        out<<endl;

        out << "var. posteriori: " << vp << endl<<endl;

        mvcXa.print(out,"MVC dos parametros ajustados");
        out<<endl;

        mvcLa.print(out,"MVC das observações ajustadas");
        out<<endl;

    }
//cout<< res<<endl;
    return res;
}

void leasq2ParallelPlanes(vec4 *plane1,vec4 *plane2,mat points1,mat points2,mat *ptOnPlane1,mat *ptOnPlane2,string repNam,uword n=5,double e=.000001)
{
    vec5 Xo,X,Xa;
    mat data,P,A,B,M,mvcXa,mvcLa,mvcV,mvcW,ptsOnPlanes;
    vec Lb,W,K,V,La,temp5,temp6;
    double vp;
    vec4 temp1,temp2,temp3,temp4;
    temp1.ones();
    temp2.zeros();
    uword iter;

    ofstream out(repNam);

    Xo.rows(0,3) = plane3points(trans(points1.row(0)),trans(points1.row(1)),trans(points1.row(2)) );
    Xo(4) = - dot(Xo.rows(0,2),points2.row(0));

    //vetor Lb

    data = join_vert(points1,points2);

    Lb = vectorise(data.t());

    //matriz dos pesos
    P = eye(Lb.n_elem,Lb.n_elem);

    //matriz jacobiana A (parametros)

    A = zeros<mat>(data.n_rows,Xo.n_elem);
    B = zeros<mat>(data.n_rows,Lb.n_elem);
    W = zeros<vec>(data.n_rows);

    A.submat(0,0,data.n_rows-1,data.n_cols-1) = data;


    for (uword i = 0; i < data.n_rows; i++)
    {
        if (i < points1.n_rows)
        {
            A(i,3) = 1;
        }
        else
        {
            A(i,4) = 1;
        }
    }

    for(uword it = 0; it < n; it ++)
    {


        uword j = 0;
        for (uword i2 = 0; i2 < data.n_rows; i2++)
        {
            B(i2,0+j)=Xo(0);
            B(i2,1+j)=Xo(1);
            B(i2,2+j)=Xo(2);

            j += 3;

            temp1.rows(0,2) = trans( data.row(i2) );
            temp2.rows(0,2) = Xo.rows(0,2);

            if(i2 < points1.n_rows)
            {
                temp2(3)=Xo(3);
            }
            else
            {
                temp2(3)=Xo(4);
            }

            W(i2) = dot(temp1,temp2);
        }

        M = B * P.i() * B.t();

        X = - inv(A.t()*M.i()*A)*(A.t()*M.i()*W);

        Xa = Xo + X;

        Xo = Xa;

        iter = it;

        if (arma::max(arma::abs(X)) < e)
        {
            break;
        }

    }

    K = -M.i()*(A*X+W);

    V = P.i() * B.t() * K;

    La = Lb+V;

    ptsOnPlanes = trans(reshape(La,data.n_cols,data.n_rows));

    *ptOnPlane1 = ptsOnPlanes.rows(0,points1.n_rows-1);

    *ptOnPlane2 = ptsOnPlanes.rows(points1.n_rows,data.n_rows-1);

    temp3.rows(0,2) =   normalise(Xa.rows(0,2));
    temp3(3)= - dot(La.rows(0,2),normalise(Xa.rows(0,2)));

    temp4 = temp3;
    temp4(3)= - dot(La.rows(La.n_elem-3,La.n_elem-1),normalise(Xa.rows(0,2)));

    *plane1 = temp3;
    *plane2 = temp4;

    vp = as_scalar( V.t() * P * V / (Lb.n_elem - Xo.n_elem) );

    mvcXa = vp * inv(A.t()*M.i()*A);

    mvcLa = vp*(inv(P)+inv(P)*B.t()*inv(M)*A*inv(A.t()*inv(M)*A)*A.t()*inv(M)*B*inv(P)-inv(P)*B.t()*inv(M)*B*inv(P));

    mvcV  = vp*P.i()-mvcLa;

    mvcW  = vp*M;

    out.precision(20);

    points1.raw_print(out,"dados de entrada (primeiro plano): ");
    out<<endl;

    points2.raw_print(out,"dados de entrada (segundo plano): ");
    out<<endl;

    out << iter <<" iterações necessárias" <<endl<<endl;

    Xa.raw_print(out," vetor Xa (pré-normalização)");
    out<<endl;

    temp3.raw_print(out,"coeficientes do primeiro plano estimado: ");
    out<<endl;

    temp4.raw_print(out,"coeficientes do segundo plano estimado: ");
    out<<endl;

    V.raw_print(out,"resíduos: ");
    out<<endl;

    ptsOnPlanes.raw_print(out,"La:");
    out<<endl;

    out << "var. posteriori: " << vp << endl<<endl;

    mvcXa.raw_print(out,"MVC dos parametros ajustados");
    out<<endl;

    mvcLa.raw_print(out,"MVC das observações ajustadas");
    out<<endl;


}

struct Sphere
{
    vec3 center;
    double r;
    mat onSpherePts;

    Sphere(mat points,string repNam,uword n=20,double e=.000001);
};

Sphere::Sphere(mat data,string repNam,uword n,double e)
{
    vec4 Xo,X,Xa;
    mat P,A,B,M,mvcXa,mvcLa,mvcV,mvcW,ptsOnPlanes;
    vec Lb,W,K,V,La,temp5,temp6,aux;
    double vp;
    vec4 temp1,temp2,temp3,temp4;
    temp1.ones();
    temp2.zeros();
    uword iter;

    ofstream out(repNam);

    aux=zeros(data.n_rows);

    for (uword i = 0; i < aux.n_rows; i++)
    {
        aux(i) = norm(data.row(i) - mean(data));
    }


    Xo.rows(1,3) = trans(mean(data));
    Xo(0) = mean(aux);

    //vetor Lb

    Lb = vectorise(data.t());

    //matriz dos pesos
    P = eye(Lb.n_elem,Lb.n_elem);

    //matriz jacobiana A (parametros)

    A = zeros<mat>(data.n_rows,Xo.n_elem);
    B = zeros<mat>(data.n_rows,Lb.n_elem);
    W = zeros<vec>(data.n_rows);

    for(uword it = 0; it < n; it ++)
    {


        uword j = 0;
        for (uword i2 = 0; i2 < data.n_rows; i2++)
        {
            A(i2,0) = - 2 * Xo(0);

            A(i2,1) = 2 * (Xo(1) - data(i2,0));
            A(i2,2) = 2 * (Xo(2) - data(i2,1));
            A(i2,3) = 2 * (Xo(3) - data(i2,2));

            B(i2,j)  = - 2 * (Xo(1) - data(i2,0));
            B(i2,1+j)= - 2 * (Xo(2) - data(i2,1));
            B(i2,2+j)= - 2 * (Xo(3) - data(i2,2));

            j += 3;


            W(i2) = std::pow((data(i2,0) - Xo(1)),2) + std::pow((data(i2,1) - Xo(2)),2) +
                    std::pow((data(i2,2) - Xo(3)),2) - std::pow(Xo(0),2);
        }

        M = B * P.i() * B.t();

        X = - inv(A.t()*M.i()*A)*(A.t()*M.i()*W);

        Xa = Xo + X;

        Xo = Xa;

        iter = it;

        if (arma::max(arma::abs(X)) < e)
        {
            break;
        }

    }

    K = -M.i()*(A*X+W);

    V = P.i() * B.t() * K;

    La = Lb+V;

    onSpherePts = trans(reshape(La,data.n_cols,data.n_rows));

    vp = as_scalar( V.t() * P * V / (Lb.n_elem - Xo.n_elem) );

    mvcXa = vp * inv(A.t()*M.i()*A);

    mvcLa = vp*(inv(P)+inv(P)*B.t()*inv(M)*A*inv(A.t()*inv(M)*A)*A.t()*inv(M)*B*inv(P)-inv(P)*B.t()*inv(M)*B*inv(P));

    mvcV  = vp*P.i()-mvcLa;

    mvcW  = vp*M;

    out.precision(20);

    data.raw_print(out,"dados de entrada :");
    out<<endl;

    out << iter <<" iterações necessárias" <<endl<<endl;

    Xa.raw_print(out,"parametros ajustados (r,XC,YC,ZC)");
    out<<endl;

    V.raw_print(out,"resíduos: ");
    out<<endl;

    La.raw_print(out,"La:");
    out<<endl;

    out << "var. posteriori: " << vp << endl<<endl;

    mvcXa.raw_print(out,"MVC dos parametros ajustados");
    out<<endl;

    mvcLa.raw_print(out,"MVC das observações ajustadas");
    out<<endl;

    //assignment to the object variables
    r = Xo(0);
    center = Xo.rows(1,3);
}

vector<vector<string>> SMMTleverARM::readLAFile(string filename)
{
//to read a file of the raw data of the Lever Arm surveying
    ifstream in(filename);
    vector<vector<string>> res;
    vector<string> A,B,C,D,E;
//to control
    int state = 0;
    string line;

    while (std::getline(in,line))
    {
        if (line.find("COD:[A") != -1)
        {
            state = 1;
            continue;
        }
        if (line.find("COD:[B") != -1)
        {
            state = 2;
            continue;
        }
        if (line.find("COD:[C") != -1)
        {
            state =3 ;
            continue;
        }
        if (line.find("COD:[D") != -1)
        {
            state =4 ;
            continue;
        }
        if (line.find("COD:[E") != -1)
        {
            state =5 ;
            continue;
        }

        if (state == 1)
        {
            A.push_back(line);
            //cout << stoi(line.substr(23,3))<<endl;
        }
        if (state == 2)
        {
            B.push_back(line);
        }
        if (state == 3)
        {
            C.push_back(line);
        }
        if (state ==4 )
        {
            D.push_back(line);
        }
        if (state == 5)
        {
            E.push_back(line);
        }
    }

    res.push_back(A);
    res.push_back(B);
    res.push_back(C);
    res.push_back(D);
    res.push_back(E);
    //cout<<A.size()<<" "<< B.size()<<" "<<C.size()<<" "<<D.size()<<" "<<E.size()<<endl;
    //cout<<res.size();

    return res;
}

SMMTleverARM::SMMTleverARM (string fileName)
{
// the constructor, all the modifying functions must have their calls here.
    data = readLAFile(fileName);
    horiz();
    zenit();
    distanc();
    avgOBS();

//the two planes
    Vplan = leasqPlane2(points.at(1),&ptsV,"plano_vert.txt",2);
    Hplan = leasqPlane2(points.at(2),&ptsH,"plano_hor.txt",2);
//the angle between them
    planesAngle = std::acos(dot(Vplan.rows(0,2),Hplan.rows(0,2))) * (180/datum::pi);

//the marks on the IMU surface
    Hm = trans(ptsH.row(2));
    Vm = trans(ptsV.row(2));

//the normal of the Hplan, is also the third direction (w),
//but at first, is needed some test
    w = Hplan.rows(0,2);
    if(w(2)<0) //the z component, needs to be positive, due to the z axis of the total station
    {
        w *= -1;
    }

//the aux1 plane
//aux1.rows(0,2) = w;
//aux1(3) = - dot(w,Vm);

//IMU center with the intersection of the plane (aux1) and the line (point Hm and direction w)
//imuC = Hm + (dot(w,(Vm-Hm)))* w; //do not delete the old implementation
    imuC = ProjPtOrtPlane(Hm,w,Vm);

//the first principal direction (u)
    u = normalise(Vm - imuC);

//the second principal direction(v)
    v = normalise(cross(w,u));

//first, translating the origin to the calculated center of the IMU
    TptMat(&ptsH,-imuC);
    TptMat(&ptsV,-imuC);
    translateAll(-imuC);
    pointsPrint("points1.txt");

//composing the first rotation matrice
    Rot1.col(0)=u;
    Rot1.col(1)=v;
    Rot1.col(2)=w;
    Rot1 = Rot1.t();

//rotating all to the imu BF
    RptMat(&ptsH,Rot1);
    RptMat(&ptsV,Rot1);
    rotAll(Rot1);
    pointsPrint("points2.txt");

//imuC += -imuC;

//creating matrices with the points, from topoPoints;
    ptsCDin  = topoPointsMat(points.at(3),0,4);
    ptsCDout = topoPointsMat(points.at(3),5,9);
    ptsCEin  = topoPointsMat(points.at(4),0,4);
    ptsCEout = topoPointsMat(points.at(4),5,9);

//the pair of parallel planes, one for each camera
    leasq2ParallelPlanes(&CDin,&CDout,ptsCDin,ptsCDout,&ptsOnCDin,&ptsOnCDout,"RcamPlanes.txt");
    leasq2ParallelPlanes(&CEin,&CEout,ptsCEin,ptsCEout,&ptsOnCEin,&ptsOnCEout,"LcamPlanes.txt");

//creating a sphere for each camera, to measure the center
    Sphere Rsphere(join_vert(ptsOnCDin,ptsOnCDout),"Rsphere.txt");
    Sphere Lsphere(join_vert(ptsOnCEin,ptsOnCEout),"Lsphere.txt");

//Y of each camera direction (with same convention)
//due to the IMU bf and the cameras position, the y component needs to be positive
// to point forward the camera axis

    yE = normalise(CEin.rows(0,2));
    if (yE(1) < 0 )
    {
        yE *= -1;
    }

    yD = normalise(CDin.rows(0,2));
    if (yD(1) < 0 )
    {
        yD *= -1;
    }

//yD.print();cout<<endl;
//yE.print();cout<<endl;

//we'll need points perfectly on the plane, and to have this guarantee:
//each one are obtained by the plane equation, with z=f(x,y)
    aux1d = onPlanePoint(CDin,ptsOnCDin(ptsOnCDin.n_rows-1,0),ptsOnCDin(ptsOnCDin.n_rows-1,1));
    aux1e = onPlanePoint(CEin,ptsOnCEin(ptsOnCEin.n_rows-1,0),ptsOnCEin(ptsOnCEin.n_rows-1,1));

    aux_D0 = onPlanePoint(CDin,ptsOnCDin(0,0),ptsOnCDin(0,1));
    aux_E0 = onPlanePoint(CEin,ptsOnCEin(0,0),ptsOnCEin(0,1));

//projecting the center of spheres on the inner plane of the rings on the camera lenses
    tempD = ProjPtOrtPlane(Rsphere.center,yD,aux1d);
    tempE = ProjPtOrtPlane(Lsphere.center,yE,aux1e);

//calculating the z direction, for each camera
    zE = normalise(aux_E0 - tempE);
    zD = normalise(aux_D0 - tempD);

//and, finally the x direction
    xE = normalise(cross(yE,zE));
    xD = normalise(cross(yD,zD));

    vec3 Xe2,Ye2,Ze2,Xd2,Yd2,Zd2;

//////formerly:
//////we can describe the axes in the camera conventional convention (z axis pointing from the CP to the focal plane)
////Xe = xE;Ye = zE;Ze = -yE;
////Xd = xD;Yd = zD;Zd = -yD;

//we can describe the axes in the camera modern convention (z axis pointing from the principle point to CP)
    Xe = xE;
    Ye = -zE;
    Ze = yE;

    Xd = xD;
    Yd = -zD;
    Zd = yD;

    //conventional convention
    Xe2 = xE;Ye2 = zE;Ze2 = -yE;
    Xd2 = xD;Yd2 = zD;Zd2 = -yD;


//with the axis, we can have the boresight matrices:
    //first with the same convention of the IMU:
    Rimu_LC2.row(0) = trans(xE);
    Rimu_LC2.row(1) = trans(yE);
    Rimu_LC2.row(2) = trans(zE);

    Rimu_RC2.row(0) = trans(xD);
    Rimu_RC2.row(1) = trans(yD);
    Rimu_RC2.row(2) = trans(zD);

    //with the camera modern convention
    Rimu_LC.row(0) = trans(Xe);
    Rimu_LC.row(1) = trans(Ye);
    Rimu_LC.row(2) = trans(Ze);

    Rimu_RC.row(0) = trans(Xd);
    Rimu_RC.row(1) = trans(Yd);
    Rimu_RC.row(2) = trans(Zd);

    //then with the camera conventional convention
    Rimu_LC3.row(0) = trans(Xe2);
    Rimu_LC3.row(1) = trans(Ye2);
    Rimu_LC3.row(2) = trans(Ze2);

    Rimu_RC3.row(0) = trans(Xd2);
    Rimu_RC3.row(1) = trans(Yd2);
    Rimu_RC3.row(2) = trans(Zd2);

//now, the lever-arm

    //antenna lever-arm: just the point at the antena top, and the offset from there to the PhC
    antLA = (points.at(0).at(0).XYZ) + PhC;

    //the lever-arm of the two cameras:
    RcamLA = Rsphere.center;
    LcamLA = Lsphere.center;

    report();
}


void SMMTleverARM::avgOBS()
{
    vector<topoPoint> temp;

    double hA,zA,dA=0,hB,zB,dB,hC,zC,dC,hD,zD,dD,hE,zE,dE,hF,zF,dF,hG,zG,dG,hH,zH,dH,hI,zI,dI,hJ,zJ,dJ;

//std::setprecision(8);
    for (unsigned int i = 0; i < hor.size(); i++)
    {
        if (i == 0)
        {
            hA = (pdpiH(hor[i][0].gdec(),hor[i][1].gdec()) +
                  pdpiH(hor[i][2].gdec(),hor[i][3].gdec()) +
                  pdpiH(hor[i][4].gdec(),hor[i][5].gdec())) / 3;

            zA = (pdpiV(zen[i][0].gdec(),zen[i][1].gdec()) +
                  pdpiV(zen[i][2].gdec(),zen[i][3].gdec()) +
                  pdpiV(zen[i][4].gdec(),zen[i][5].gdec())) / 3;

            for (unsigned int j = 0; j < dists[i].size(); j++)
            {
                dA += dists[i][j]/6;
            }

            topoPoint ANT(hA,zA,dA,false);
            temp.push_back(ANT);

            //cout <<hA<<" "<<zA<<" "<<dA<<endl;
            //cout <<hA<<" "<<zA<<" "<<dA<<endl;
        }
        if (i == 1)
        {
            hA = (pdpiH(hor[i][0].gdec(),hor[i][5].gdec()) +
                  pdpiH(hor[i][10].gdec(),hor[i][15].gdec())) / 2;

            hB = (pdpiH(hor[i][1].gdec(),hor[i][6].gdec()) +
                  pdpiH(hor[i][11].gdec(),hor[i][16].gdec())) / 2;

            hC = (pdpiH(hor[i][2].gdec(),hor[i][7].gdec()) +
                  pdpiH(hor[i][12].gdec(),hor[i][17].gdec())) / 2;

            hD = (pdpiH(hor[i][3].gdec(),hor[i][8].gdec()) +
                  pdpiH(hor[i][13].gdec(),hor[i][18].gdec())) / 2;

            hE = (pdpiH(hor[i][4].gdec(),hor[i][9].gdec()) +
                  pdpiH(hor[i][14].gdec(),hor[i][19].gdec())) / 2;



            zA = (pdpiV(zen[i][0].gdec(),zen[i][5].gdec()) +
                  pdpiV(zen[i][10].gdec(),zen[i][15].gdec())) / 2;

            zB = (pdpiV(zen[i][1].gdec(),zen[i][6].gdec()) +
                  pdpiV(zen[i][11].gdec(),zen[i][16].gdec())) / 2;

            zC = (pdpiV(zen[i][2].gdec(),zen[i][7].gdec()) +
                  pdpiV(zen[i][12].gdec(),zen[i][17].gdec())) / 2;

            zD = (pdpiV(zen[i][3].gdec(),zen[i][8].gdec()) +
                  pdpiV(zen[i][13].gdec(),zen[i][18].gdec())) / 2;

            zE = (pdpiV(zen[i][4].gdec(),zen[i][9].gdec()) +
                  pdpiV(zen[i][14].gdec(),zen[i][19].gdec())) / 2;


            dA = (dists[i][0]+dists[i][5]+dists[i][10]+dists[i][15]) / 4;

            dB = (dists[i][1]+dists[i][6]+dists[i][11]+dists[i][16]) / 4;

            dC = (dists[i][2]+dists[i][7]+dists[i][12]+dists[i][17]) / 4;

            dD = (dists[i][3]+dists[i][8]+dists[i][13]+dists[i][18]) / 4;

            dE = (dists[i][4]+dists[i][9]+dists[i][14]+dists[i][19]) / 4;


            topoPoint A(hA,zA,dA,false);
            temp.push_back(A);
            topoPoint B(hB,zB,dB,false);
            temp.push_back(B);
            topoPoint C(hC,zC,dC,false);
            temp.push_back(C);
            topoPoint D(hD,zD,dD,false);
            temp.push_back(D);
            topoPoint E(hE,zE,dE,false);
            temp.push_back(E);

//                                cout <<hA<<" "<<zA<<" "<<dA<<endl;


        }
        if (i == 2)
        {
            hA = (pdpiH(hor[i][0].gdec(),hor[i][5].gdec()) +
                  pdpiH(hor[i][10].gdec(),hor[i][15].gdec())) / 2;

            hB = (pdpiH(hor[i][1].gdec(),hor[i][6].gdec()) +
                  pdpiH(hor[i][11].gdec(),hor[i][16].gdec())) / 2;

            hC = (pdpiH(hor[i][2].gdec(),hor[i][7].gdec()) +
                  pdpiH(hor[i][12].gdec(),hor[i][17].gdec())) / 2;

            hD = (pdpiH(hor[i][3].gdec(),hor[i][8].gdec()) +
                  pdpiH(hor[i][13].gdec(),hor[i][18].gdec())) / 2;

            hE = (pdpiH(hor[i][4].gdec(),hor[i][9].gdec()) +
                  pdpiH(hor[i][14].gdec(),hor[i][19].gdec())) / 2;



            zA = (pdpiV(zen[i][0].gdec(),zen[i][5].gdec()) +
                  pdpiV(zen[i][10].gdec(),zen[i][15].gdec())) / 2;

            zB = (pdpiV(zen[i][1].gdec(),zen[i][6].gdec()) +
                  pdpiV(zen[i][11].gdec(),zen[i][16].gdec())) / 2;

            zC = (pdpiV(zen[i][2].gdec(),zen[i][7].gdec()) +
                  pdpiV(zen[i][12].gdec(),zen[i][17].gdec())) / 2;

            zD = (pdpiV(zen[i][3].gdec(),zen[i][8].gdec()) +
                  pdpiV(zen[i][13].gdec(),zen[i][18].gdec())) / 2;

            zE = (pdpiV(zen[i][4].gdec(),zen[i][9].gdec()) +
                  pdpiV(zen[i][14].gdec(),zen[i][19].gdec())) / 2;


            dA = (dists[i][0]+dists[i][5]+dists[i][10]+dists[i][15]) / 4;

            dB = (dists[i][1]+dists[i][6]+dists[i][11]+dists[i][16]) / 4;

            dC = (dists[i][2]+dists[i][7]+dists[i][12]+dists[i][17]) / 4;

            dD = (dists[i][3]+dists[i][8]+dists[i][13]+dists[i][18]) / 4;

            dE = (dists[i][4]+dists[i][9]+dists[i][14]+dists[i][19]) / 4;


            topoPoint A(hA,zA,dA,false);
            temp.push_back(A);
            topoPoint B(hB,zB,dB,false);
            temp.push_back(B);
            topoPoint C(hC,zC,dC,false);
            temp.push_back(C);
            topoPoint D(hD,zD,dD,false);
            temp.push_back(D);
            topoPoint E(hE,zE,dE,false);
            temp.push_back(E);

//                                cout <<hA<<" "<<zA<<" "<<dA<<endl;

        }
        if (i == 3)
        {
            topoPoint EXTRA(hor[i][0].gdec(),zen[i][0].gdec(),dists[i][0],false);

            hA = pdpiH(hor[i][1].gdec(),hor[i][11].gdec());
            hB = pdpiH(hor[i][2].gdec(),hor[i][12].gdec());
            hC = pdpiH(hor[i][3].gdec(),hor[i][13].gdec());
            hD = pdpiH(hor[i][4].gdec(),hor[i][14].gdec());
            hE = pdpiH(hor[i][5].gdec(),hor[i][15].gdec());
            hF = pdpiH(hor[i][6].gdec(),hor[i][16].gdec());
            hG = pdpiH(hor[i][7].gdec(),hor[i][17].gdec());
            hH = pdpiH(hor[i][8].gdec(),hor[i][18].gdec());
            hI = pdpiH(hor[i][9].gdec(),hor[i][19].gdec());
            hJ = pdpiH(hor[i][10].gdec(),hor[i][20].gdec());

            zA = pdpiV(zen[i][1].gdec(),zen[i][11].gdec());
            zB = pdpiV(zen[i][2].gdec(),zen[i][12].gdec());
            zC = pdpiV(zen[i][3].gdec(),zen[i][13].gdec());
            zD = pdpiV(zen[i][4].gdec(),zen[i][14].gdec());
            zE = pdpiV(zen[i][5].gdec(),zen[i][15].gdec());
            zF = pdpiV(zen[i][6].gdec(),zen[i][16].gdec());
            zG = pdpiV(zen[i][7].gdec(),zen[i][17].gdec());
            zH = pdpiV(zen[i][8].gdec(),zen[i][18].gdec());
            zI = pdpiV(zen[i][9].gdec(),zen[i][19].gdec());
            zJ = pdpiV(zen[i][10].gdec(),zen[i][20].gdec());

            dA = (dists[i][1]+dists[i][11]) / 2;
            dB = (dists[i][2]+dists[i][12]) / 2;
            dC = (dists[i][3]+dists[i][13]) / 2;
            dD = (dists[i][4]+dists[i][14]) / 2;
            dE = (dists[i][5]+dists[i][15]) / 2;
            dF = (dists[i][6]+dists[i][16]) / 2;
            dG = (dists[i][7]+dists[i][17]) / 2;
            dH = (dists[i][8]+dists[i][18]) / 2;
            dI = (dists[i][9]+dists[i][19]) / 2;
            dJ = (dists[i][10]+dists[i][20]) / 2;

            topoPoint A(hA,zA,dA,false);
            temp.push_back(A);
            topoPoint B(hB,zB,dB,false);
            temp.push_back(B);
            topoPoint C(hC,zC,dC,false);
            temp.push_back(C);
            topoPoint D(hD,zD,dD,false);
            temp.push_back(D);
            topoPoint E(hE,zE,dE,false);
            temp.push_back(E);
            topoPoint F(hF,zF,dF,false);
            temp.push_back(F);
            topoPoint G(hG,zG,dG,false);
            temp.push_back(G);
            topoPoint H(hH,zH,dH,false);
            temp.push_back(H);
            topoPoint I(hI,zI,dI,false);
            temp.push_back(I);
            topoPoint J(hJ,zJ,dJ,false);
            temp.push_back(J);

            temp.push_back(EXTRA);

//                               cout <<hA<<" "<<zA<<" "<<dA<<endl;


        }
        if (i == 4)
        {
            topoPoint EXTRA(hor[i][0].gdec(),zen[i][0].gdec(),dists[i][0],false);

            hA = pdpiH(hor[i][1].gdec(),hor[i][11].gdec());
            hB = pdpiH(hor[i][2].gdec(),hor[i][12].gdec());
            hC = pdpiH(hor[i][3].gdec(),hor[i][13].gdec());
            hD = pdpiH(hor[i][4].gdec(),hor[i][14].gdec());
            hE = pdpiH(hor[i][5].gdec(),hor[i][15].gdec());
            hF = pdpiH(hor[i][6].gdec(),hor[i][16].gdec());
            hG = pdpiH(hor[i][7].gdec(),hor[i][17].gdec());
            hH = pdpiH(hor[i][8].gdec(),hor[i][18].gdec());
            hI = pdpiH(hor[i][9].gdec(),hor[i][19].gdec());
            hJ = pdpiH(hor[i][10].gdec(),hor[i][20].gdec());

            zA = pdpiV(zen[i][1].gdec(),zen[i][11].gdec());
            zB = pdpiV(zen[i][2].gdec(),zen[i][12].gdec());
            zC = pdpiV(zen[i][3].gdec(),zen[i][13].gdec());
            zD = pdpiV(zen[i][4].gdec(),zen[i][14].gdec());
            zE = pdpiV(zen[i][5].gdec(),zen[i][15].gdec());
            zF = pdpiV(zen[i][6].gdec(),zen[i][16].gdec());
            zG = pdpiV(zen[i][7].gdec(),zen[i][17].gdec());
            zH = pdpiV(zen[i][8].gdec(),zen[i][18].gdec());
            zI = pdpiV(zen[i][9].gdec(),zen[i][19].gdec());
            zJ = pdpiV(zen[i][10].gdec(),zen[i][20].gdec());

            dA = (dists[i][1]+dists[i][11]) / 2;
            dB = (dists[i][2]+dists[i][12]) / 2;
            dC = (dists[i][3]+dists[i][13]) / 2;
            dD = (dists[i][4]+dists[i][14]) / 2;
            dE = (dists[i][5]+dists[i][15]) / 2;
            dF = (dists[i][6]+dists[i][16]) / 2;
            dG = (dists[i][7]+dists[i][17]) / 2;
            dH = (dists[i][8]+dists[i][18]) / 2;
            dI = (dists[i][9]+dists[i][19]) / 2;
            dJ = (dists[i][10]+dists[i][20]) / 2;

            topoPoint A(hA,zA,dA,false);
            temp.push_back(A);
            topoPoint B(hB,zB,dB,false);
            temp.push_back(B);
            topoPoint C(hC,zC,dC,false);
            temp.push_back(C);
            topoPoint D(hD,zD,dD,false);
            temp.push_back(D);
            topoPoint E(hE,zE,dE,false);
            temp.push_back(E);
            topoPoint F(hF,zF,dF,false);
            temp.push_back(F);
            topoPoint G(hG,zG,dG,false);
            temp.push_back(G);
            topoPoint H(hH,zH,dH,false);
            temp.push_back(H);
            topoPoint I(hI,zI,dI,false);
            temp.push_back(I);
            topoPoint J(hJ,zJ,dJ,false);
            temp.push_back(J);

            temp.push_back(EXTRA);

//                               cout <<hA<<" "<<zA<<" "<<dA<<endl;

        }


        points.push_back(temp);
        temp.clear();
    }

    pointsPrint("points0.txt");
}

void SMMTleverARM::pointsPrint(string filename)
{
    ofstream outFile(filename);

    for (unsigned ii = 0; ii < points.size(); ii++)
    {
        for (unsigned jj = 0; jj < points[ii].size(); jj++)
        {
            outFile << points[ii][jj].XYZ(0) <<" ";
            outFile << points[ii][jj].XYZ(1) <<" ";
            outFile << points[ii][jj].XYZ(2) << endl;
        }
    }

    outFile.close();
}


void SMMTleverARM::horiz()
{
    vector<GMS> temp;
    int g,m,s;

    for (unsigned int i = 0; i < data.size() ; i++)
    {

        for (unsigned int j = 0; j < data[i].size(); j++)
        {

            g = stoi(data[i][j].substr(23,3));
            m = stoi(data[i][j].substr(27,2));
            s = stoi(data[i][j].substr(29,2));
            GMS temp2(g,m,s);

            temp.push_back(temp2);
        }
        hor.push_back(temp);
        temp.clear();
    }

}

void SMMTleverARM::zenit()
{
    vector<GMS> temp;
    int g,m,s;

    for (unsigned int i = 0; i < data.size() ; i++)
    {

        for (unsigned int j = 0; j < data[i].size(); j++)
        {

            g = stoi(data[i][j].substr(35,3));
            m = stoi(data[i][j].substr(39,2));
            s = stoi(data[i][j].substr(41,2));
            GMS temp2(g,m,s);
//            cout<<temp2.gdec()<<endl;

            temp.push_back(temp2);
        }
        zen.push_back(temp);
        temp.clear();
    }

}

void SMMTleverARM::distanc()
{
    vector<double> temp;
    double temp2;

    for (unsigned int i = 0; i < data.size() ; i++)
    {

        for (unsigned int j = 0; j < data[i].size(); j++)
        {

            temp2 = stod(data[i][j].substr(47,10));
//            cout<<temp2<<endl;

            temp.push_back(temp2);
        }
        dists.push_back(temp);
        temp.clear();
    }

}

topoPoint::topoPoint(double h,double z,double di,bool isRad)
{
    Hz=h;
    Ze=z;
    Di=di;
    rad = isRad;
    cart();
}

void topoPoint::cart()
{
    double X,Y,Z;

    if (rad)
    {
        X = std::sin(Ze) * std::sin(Hz) * Di;
        Y = std::sin(Ze) * std::cos(Hz) * Di;
        Z = std::cos(Ze) * Di;
    }
    else
    {
        X = std::sin(Ze*(datum::pi/180)) * std::sin(Hz*(datum::pi/180)) * Di;
        Y = std::sin(Ze*(datum::pi/180)) * std::cos(Hz*(datum::pi/180)) * Di;
        Z = std::cos(Ze*(datum::pi/180)) * Di;
    }

    vec3 temp = {X,Y,Z};

    XYZ = temp;
}

double pdpiH(double h1,double h2)
{
    if (h1 > h2)
    {
        return (h1+h2+180)/2;
    }
    else
    {
        return (h1+h2-180)/2;
    }
}

double pdpiV(double v1,double v2)
{
    if (v1 < v2)
    {
        return (v1-v2+360)/2;
    }
    else
    {
        return (v2-v1+360)/2;
    }
}

vec4 plane3points(vec3 u,vec3 v,vec3 w)
{
//the function returns a 4-vector with the plane equation coefficients
//assuming the form  ax + by + cz + d = 0
    vec4 res;
    vec3 p1,p2,n;
    vec d;
    p1 = normalise(v - u);
    p2 = normalise(w - u);
    n  = normalise(cross(p1,p2));
    d = - dot(n,u);

    res = join_vert(n,d);

//cout << res <<endl;
    return res;
}

void topoPoint::translate(vec3 tvec)
{
    XYZ += tvec;
}

void topoPoint::Rotate(mat R)
{
    XYZ = R * XYZ;
}

void TptMat (mat *M,vec Tvec)
{
    mat temp = *M;

    if(temp.n_cols == Tvec.n_elem)
    {
        for (uword i = 0; i < temp.n_rows; i++)
        {
            temp.row(i) = trans(trans(temp.row(i)) + Tvec );
        }
    }
    else
    {
        cout << "nada realizado, tamanho de colunas da matriz diferente do numero de elementos do vetor dado"<<endl;
    }

    *M = temp;
}

void RptMat (mat *M,mat R)
{
    mat temp = *M;

    if( (temp.n_cols == R.n_cols) & (R.n_cols == R.n_rows) )
    {
        for (uword i = 0; i < temp.n_rows; i++)
        {
            temp.row(i) = trans(R*trans(temp.row(i)));
        }
    }
    else
    {
        cout << "nada realizado, tamanho de colunas da matriz diferente do numero de colunas da matriz dada"<<endl;
        cout<<"ou a matriz fornecida não é quadrada"<<endl;
    }

    *M = temp;
}

void SMMTleverARM::translateAll(vec3 tvec)
{
    for (unsigned int i = 0; i < points.size(); i++)
    {
        for (unsigned int j = 0; j < points[i].size(); j++)
        {
            points[i][j].translate(tvec);
        }
    }
}

void SMMTleverARM::rotAll (mat R)
{
    for (unsigned int i = 0; i < points.size(); i++)
    {
        for (unsigned int j = 0; j < points[i].size(); j++)
        {
            points[i][j].Rotate(R);
        }
    }
}

mat topoPointsMat (vector<topoPoint> points,uword first,uword last)
{
    mat res;

    uword nlin = last - first + 1;

    res.zeros(nlin,3);

    for (uword i = first; i < last+1; i++)
    {
        res.row(i-first) = trans(points.at(i).XYZ);
//        cout << i <<endl;
    }

//    cout<<res<<endl;

    return res;
}

vec4 parallelPlane(vec4 plane1,vec3 pointOnPlane)
{
    vec4 res;
    vec3 n = normalise(plane1.rows(0,2));

    res.rows(0,2) = n;
    res(3) = - dot(n,pointOnPlane);

    return res;
}

vec3 ProjPtOrtPlane(vec3 pt2Proj,vec3 plNorm,vec3 ptOnPlane)
{
    return pt2Proj + (dot(normalise(plNorm),ptOnPlane-pt2Proj)) * normalise(plNorm);
}

vec3 onPlanePoint(vec4 plCoef,double x,double y)
{
    vec3 res;

    res(0) = x;
    res(1) = y;
    res(2) = - (plCoef(0)*x+plCoef(1)*y+plCoef(3))/plCoef(2);

    return res;
}

void SMMTleverARM::report()
{
    ofstream out("report.txt");
    out.precision(8);

    out << "Relatório de Saída do Processamento da ";
    out<< "determinação dos parâmetros de orientação relativa"   <<endl;
    out<<"do Sistema de Mapeamento Móvel Terrestre do LAPE"<<endl<<endl<<endl;

    out << "Lever Arm (no body-frame da IMU) --------------------------------"<<endl<<endl;

    out<< "IMU -> Antena"<<endl;
    out << "X(m): " <<antLA(0)<< " | Y(m): " <<antLA(1)<< " | Z(m): " <<antLA(2)<<endl<<endl;

    out<< "IMU -> Câmera Esquerda"<<endl;
    out << "X(m): " <<LcamLA(0)<< " | Y(m): " <<LcamLA(1)<< " | Z(m): " <<LcamLA(2)<<endl<<endl;

    out<< "IMU -> Câmera Direita"<<endl;
    out << "X(m): " <<RcamLA(0)<< " | Y(m): " <<RcamLA(1)<< " | Z(m): " <<RcamLA(2)<<endl<<endl<<endl;

    out<<"Distancias:"<<endl;
    out<<"Entre câmeras: "<<arma::norm(LcamLA-RcamLA)<<endl;
    out<<"IMU - Antena: "<<arma::norm(antLA)<<endl;
    out<<"IMU - Cam Esq. : "<<arma::norm(LcamLA)<<endl;
    out<<"IMU - Cam Dir. : "<<arma::norm(RcamLA)<<endl;


    out<< "Matrizes do Boresight (cossenos diretores): "<<endl<<endl;
    Rimu_LC.raw_print(out,"IMU (BF) -> Câmera Esquerda (BF)");
    out<<endl;
    Rimu_RC.raw_print(out,"IMU (BF) -> Câmera Direita (BF)");
    out<<endl<<endl;

    out<< "Matrizes do Boresight (cossenos diretores, mesma convenção da IMU): "<<endl<<endl;
    Rimu_LC2.raw_print(out,"IMU (BF) -> Câmera Esquerda (BF)");
    out<<endl;
    Rimu_RC2.raw_print(out,"IMU (BF) -> Câmera Direita (BF)");
    out<<endl<<endl;

    out<<"no último caso, os angulos formados com os eixos do BF da IMU: "<<endl<<endl;

    out<<"IMU (BF) -> Câmera Esquerda (BF)"<<endl<<arma::acos(Rimu_LC2)*(180/datum::pi)<<endl<<endl;
    out<<"IMU (BF) -> Câmera Direita (BF)"<<endl<<arma::acos(Rimu_RC2)*(180/datum::pi)<<endl;

    out<<endl<<"testes de consistência dos resultados"<<endl;
    out<<"as matrizes de rotação deverão, ao ser multiplicadas por sua transposta, serem iguais a matriz identidade"<<endl<<endl;

    out<<"Rtopo_imu * Rtopo_imu.t()"<<endl<<Rot1 * Rot1.t()<<endl;

    out<<"Rimu_LC * Rimu_LC.t()"<<endl<<Rimu_LC * Rimu_LC.t()<<endl;
    out<<"Rimu_RC * Rimu_RC.t()"<<endl<<Rimu_RC * Rimu_RC.t()<<endl;

    out<<"Rimu_LC2 * Rimu_LC2.t()"<<endl<<Rimu_LC2 * Rimu_LC2.t()<<endl;
    out<<"Rimu_RC2 * Rimu_RC2.t()"<<endl<<Rimu_RC2 * Rimu_RC2.t()<<endl;

}

