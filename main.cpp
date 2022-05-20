#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <eigen3/Eigen/Dense>

#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;


double f1 (double t, VectorXd y) {
    double f = 2 * t * y(0) * log(max(y(1), 0.001));
    return f;
}

double f2 (double t, VectorXd y) {
    double f = -2 * t * y(1) * log(max(y(0), 0.001));
    return f;
}

double AutoStep (double h0, double tol, double err, double facmax) {
    double s = 4.;
    
    double facmin = 0.2;
    double fac = 0.8;
    
    double h = h0 * min(facmax, max(facmin, fac * pow(tol/err, 1/4) ));
    
    return h;
}

vector<VectorXd> k_coeff (double t, VectorXd y, double h) {
    vector<vector<double>> a = {{}, {1/2}, {0, 1/2}, {0, 0, 1}};
    vector<double> c = {0, 1/2, 1/2, 1};
    vector<VectorXd> k(4);
    VectorXd kv1(2);
    kv1(0) = f1(t, y);
    kv1(1) = f2(t, y);
    k[0] = kv1;
    
    VectorXd kv2(2);
    kv2(0) = f1(t + c[1] * h, y + h * a[1][0] * k[0]);
    kv2(1) = f2(t + c[1] * h, y + h * a[1][0] * k[0]);
    k[1] = kv2;
    
    VectorXd kv3(2);
    kv3(0) = f1(t + c[2] * h, h * (a[2][0] * k[0] + a[2][1] * k[1]));
    kv3(1) = f2(t + c[2] * h, h * (a[2][0] * k[0] + a[2][1] * k[1]));
    k[2] = kv3;
    
    VectorXd kv4(2);
    kv4(0) = f1(t + c[3] * h, y + h * (a[3][0] * k[0] + a[3][1] * k[1] + a[3][2] * k[2]));
    kv4(1) = f2(t + c[3] * h, y + h * (a[3][0] * k[0] + a[3][1] * k[1] + a[3][2] * k[2]));
    k[3] = kv4;
    
    return k;
}

VectorXd RK_Step (vector<VectorXd> k, VectorXd y, double h) {
    vector<double> b = {1/6, 2/6, 2/6, 1/6};
    VectorXd f = y + h * (b[0] * k[0] + b[1] * k[1] + b[2] * k[2] + b[3] * k[3]);
//    cout << (b[0] * k[0] + b[1] * k[1] + b[2] * k[2] + b[3] * k[3]) << endl;
//    cout << k[0] << " " << k[1] << " " << k[2] << " " << k[3] << endl;
//    cout << k[0] << endl;
//    cout << endl;
//    cout << k[1] << endl;
//    cout << endl;
//    cout << k[2] << endl;
//    cout << endl;
//    cout << k[3] << endl;
//    cout << endl;
    return f;
}

void RungeKutta (vector<VectorXd>& y, vector<double>& t) {
    double Nv = 6.0;
    
    
    double tol = 0.001;
    
    double t1 = Nv * 0.1;
    double t2 = t1 + 4.0;

    t.push_back(t1);
    
    double y10 = exp(sin(t1 * t1));
    double y20 = exp(cos(t1 * t1));
    VectorXd y0 {{y10, y20}};

    y.push_back(y0);
    
    double t_curr = t1;
    double h_curr = 1;
    double facmax = 1.5;
    
    while (t_curr <= t2) {
//        cout << h_curr << endl;
        while (true) {
            VectorXd y1 = RK_Step(k_coeff(t_curr, y[y.size()-1], h_curr), y[y.size()-1], h_curr);
            VectorXd y2 = RK_Step(k_coeff(t_curr + h_curr, y1, h_curr), y1, h_curr);
            VectorXd w = RK_Step(k_coeff(t_curr, y[y.size()-1], 2 * h_curr), y[y.size()-1], 2 * h_curr);
            
            
//            cout << y1 << " " << y[y.size()-1] << endl;
            
            
            double err =1 / (pow(2, 4)- 1) * (y2-w).norm();
            double h_new = AutoStep(h_curr, tol, err, facmax);
            if (err <= tol) {
                y.push_back(y2);
                t_curr += 2 * h_curr;
                t.push_back(t_curr);
                h_curr = h_new;
                break;
            } else {
                h_curr = h_new;
                facmax = 1;
            }
        }
    }
}

vector<VectorXd> k_coeffD (double t, VectorXd y, double h) {
    vector<vector<double>> a = {{}, {1/5}, {3/40, 9/40}, {44/45, -56/15, 32/9}, {19372/6561, -25360/2187, 64448/6561, -212/729}, {9017/3168, -355/33, 46732/5247, 49/176, -5103/18656}, {35/384, 0, 500/1113, 125/192, -2187/6784, 1184}};
    vector<double> c = {0, 1/5, 3/10, 4/5, 8/9, 1, 1};
    vector<VectorXd> k(4);
    VectorXd kv1(2);
    kv1(0) = f1(t, y);
    kv1(1) = f2(t, y);
    k[0] = kv1;
    
    VectorXd kv2(2);
    kv2(0) = f1(t + c[1] * h, y + h * a[1][0] * k[0]);
    kv2(1) = f2(t + c[1] * h, y + h * a[1][0] * k[0]);
    k[1] = kv2;
    
    VectorXd kv3(2);
    kv3(0) = f1(t + c[2] * h, h * (a[2][0] * k[0] + a[2][1] * k[1]));
    kv3(1) = f2(t + c[2] * h, h * (a[2][0] * k[0] + a[2][1] * k[1]));
    k[2] = kv3;
    
    VectorXd kv4(2);
    kv4(0) = f1(t + c[3] * h, y + h * (a[3][0] * k[0] + a[3][1] * k[1] + a[3][2] * k[2]));
    kv4(1) = f2(t + c[3] * h, y + h * (a[3][0] * k[0] + a[3][1] * k[1] + a[3][2] * k[2]));
    k[3] = kv4;
    
    return k;
}

VectorXd DP_Step1 (vector<VectorXd> k, VectorXd y, double h) {
    vector<double> b = {35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0};
    VectorXd y1 = y + h * (b[0] * k[0] + b[1] * k[1] + b[2] * k[2] + b[3] * k[3] + b[4] * k[4] + b[5] * k[5] + b[6] * k[6]);
    return y1;
}

VectorXd DP_Step2 (vector<VectorXd> k, VectorXd y, double h) {
    vector<double> b = {5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40};
    VectorXd y1 = y + h * (b[0] * k[0] + b[1] * k[1] + b[2] * k[2] + b[3] * k[3] + b[4] * k[4] + b[5] * k[5] + b[6] * k[6]);
    return y1;
}

void DormPrince (vector<double> k) {
    double Nv = 6.0;


    double tol = 0.001;

    double t1 = Nv * 0.1;
    double t2 = t1 + 4.0;
    vector<double> t;
    t.push_back(t1);

    double y10 = exp(sin(t1 * t1));
    double y20 = exp(cos(t1 * t1));
    VectorXd y0 {{y10, y20}};
    vector<VectorXd> y;
    y.push_back(y0);

    double t_curr = t1;
    double h_curr = 0.1;
    double facmax = 1.5;

    while (t_curr <= t2) {
        while (true) {
            VectorXd y1 = DP_Step1(k_coeffD(t_curr, *y.end(), h_curr), *y.end(), h_curr);
            VectorXd y2 = DP_Step2(k_coeffD(t_curr, *y.end(), h_curr), *y.end(), h_curr);
            double err = (y2-y1).norm();
            double h_new = AutoStep(h_curr, tol, err, facmax);
            if (err <= tol) {
                y.push_back(y1);
                t_curr += h_curr;
                t.push_back(t_curr);
                h_curr = h_new;
                break;
            } else {
                h_curr = h_new;
            }
        }
    }
}

int main(int argc, const char * argv[]) {
    vector<VectorXd> y;
    vector<double> t;
    RungeKutta(y, t);

    double max_diff = 0;
    
    for (int i = 0; i < y.size(); ++i) {
        double ti = t[i];
        VectorXd yext {{exp(sin(ti * ti)), exp(cos(ti * ti))}};
        double nor = (y[i]-yext).norm();
        if (nor > max_diff) {
            max_diff = nor;
        }
//        cout << nor << endl;
    }
    cout << max_diff << endl;

    return 0;
}