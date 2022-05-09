#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

vector<double> SturmLiuvill_Problem (vector<double> THETA, vector<double> PHI, double r, double R_k, double L, int t_max) 
{
    double theta_k00 = 1/(R_k * sqrt(M_PI * L));
    double theta_k0q = 0;
    double theta_kp0 = 0;
    double theta_kpq = 0;
    
    double  j_1q = 1, z = 1;
 
    for (int q = 0; q < Infinity; ++q) 
    {
        theta_k0q += j0(r/(R_k * j_1q))/j0(j_1q);
        for (int p = 1; p < Infinity; ++p) 
        {
            theta_kp0 += cos(p * M_PI * z / L) / (R_k * sqrt(M_PI * L / 2));
            theta_kpq += theta_kp0 * theta_k0q;
        }
    }
    for (int t = 0; t < t_max; ++t) 
    {
        THETA.push_back(theta_k00 + theta_k0q * PHI[t] + theta_kp0 * PHI[t] + theta_kpq * PHI[t]);
    }
    
    return THETA;
}

vector<double> FurieCoeff (double alpha, double beta, double a, double b, double c, double d) 
{
    vector<double> I;
    
    return I;
}

vector<double> FurieInitial (double R, double q, double p) 
{
    double J = (R * R * j1(q)) / (R * j1(q) - R * j1(p));
    
    vector<double> Ck;
    return Ck;
}

vector<double> CoordinateFunction (vector<double> PHI, int a_max, int b_max) 
{
    for (int a = 0; a < a_max; ++a) {
        for (int b = 0; b < b_max; ++b) {
            
        }
    }
    return PHI;
}

int main()
{
    vector<string> msg {"Hello", "C++", "World", "from", "VS Code", "and the C++ extension!"};

    for (const string& word : msg)
    {
        cout << word << " ";
    }
    cout << endl;
}