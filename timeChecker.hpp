#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <utility>

using namespace std;

//struct valAndPos
//{
//    int pos;
//    double val;
//
//    valAndPos(double v,int p);
//};
//
//valAndPos::valAndPos(double v,int p)
//{
//    pos = p;val = v;
//}

/// vector<pair<double,int>>

typedef vector<pair<double,int>> val_pos;
typedef vector<pair<unsigned int,unsigned int>> ipair;

val_pos listCreator(vector<double> input)
{
    val_pos res;

    for (unsigned int i = 0; i<input.size(); i++)
    {
        pair<double,int> temp(input.at(i),i);
        res.push_back(temp);
    }

    return res;
}

void checkAndRemove (val_pos &input,double val,vector<double> difs,bool minor = true)
{
    for (unsigned int i = 1; i<input.size(); i++)
    {
        if (minor)
        {
            if(difs.at(i-1) < val)
            {
                input.erase(input.begin()+i);
            }
        }
        else
        {
            if(difs.at(i-1) > val)
            {
                input.erase(input.begin()+i);
            }
        }
    }
}

vector<double> dif_vector2(val_pos a,val_pos b)
{
    vector<double> res;

    double temp;

    if (a.size() > b.size())
    {
        for (unsigned int i = 0; i<b.size(); i++)
        {
            temp = a.at(i).first - b.at(i).first;
            res.push_back(temp);
        }
    }
    else
    {
        for (unsigned int i = 0; i<a.size(); i++)
        {
            temp = a.at(i).first - b.at(i).first;
            res.push_back(temp);
        }
    }

    return res;
}


vector<double> dif_vector(vector<double> a,vector<double> b)
{
    vector<double> res;

    double temp;

    if (a.size() > b.size())
    {
        for (unsigned int i = 0; i<b.size(); i++)
        {
            temp = a.at(i) - b.at(i);
            res.push_back(temp);
        }
    }
    else
    {
        for (unsigned int i = 0; i<a.size(); i++)
        {
            temp = a.at(i) - b.at(i);
            res.push_back(temp);
        }
    }

    return res;
}

vector<int> check_invalid(vector<double> input,double val,bool minor = true)
{
//return the invalid indexes

    vector<int> res;

    for (unsigned int i = 0; i<input.size(); i++)
    {
        if (minor)
        {
            if(input.at(i)<val)
            {
                res.push_back(i);
            }
        }
        else
        {
            if(input.at(i)>val)
            {
                res.push_back(i);
            }
        }
    }

    return res;
}

void deletePositions(vector<int> positions,vector<double> &inputVector)
{
    //the indices must be ordered
    std::sort (positions.begin(),positions.end());

    int excluded = 0;
    for (unsigned int j = 0; j<positions.size(); j++)
    {
        inputVector.erase(inputVector.begin()+positions.at(j)-excluded);
        excluded++;

//                cout<<j<<" "<<positions.at(j)<<" "<<excluded<<" "<<positions.at(j-excluded)<<endl;
    }
}

struct timeFiltering
{
    vector<int> excludeA,excludeB;

    vector<double> inDifA,inDifB;

    val_pos A,B;

    bool finish = false;

    ipair finalpairs;

    timeFiltering(vector<double> a,vector<double> b);
};

bool checkDiffs(vector<double> difVector,val_pos &a,val_pos &b,double checkVal = .1)
{

    for (unsigned int i = 0; i<difVector.size(); i++)
    {
        if(fabs(difVector.at(i)) > checkVal)
        {
            if (difVector.at(i) < 0)
            {
                a.erase(a.begin()+i+1);
                return false;
            }
            else
            {
                b.erase(b.begin()+i+1);
                return false;
            }
        }
    }
    return true;
}

vector<double> internalDifs(val_pos input)
{
    vector<double> res;

    //yeah, this one begins at index number 1
    for (unsigned int i = 1; i<input.size(); i++)
    {
        double temp = input.at(i).first - input.at(i-1).first;
//        cout<<temp<<" ";
        res.push_back(temp);
    }
//    cout<<endl;

    return res;
}


timeFiltering::timeFiltering(vector<double> a,vector<double> b)
{
    int d = 0;


    A = listCreator(a);
    B = listCreator(b);


    //calculating the internal diferences
    inDifA = internalDifs(A);
    inDifB = internalDifs(B);

        cout<<A.size()<<"  ";
        cout<<B.size()<<endl;

    //first: checking if there's invalid values in the input vectors
    checkAndRemove(A,0.2,inDifA);
    checkAndRemove(B,0.2,inDifB);

        cout<<A.size()<<"  ";
        cout<<B.size()<<endl<<endl;


    for(;;)
    {


        //calculating again
        inDifA = internalDifs(A);
        inDifB = internalDifs(B);



        vector<double> diffs = dif_vector(inDifA,inDifB);



        finish = checkDiffs(diffs,A,B);



        if (finish && A.size() == B.size())
        {
            break;
        }

        cout<<A.size()<<"  ";
        cout<<B.size()<<endl;

    }

    for (unsigned int i = 0;i<A.size();i++)
    {

        pair<unsigned int,unsigned int> temp(A.at(i).second,B.at(i).second);
        finalpairs.push_back(temp);

    }

    //
}
