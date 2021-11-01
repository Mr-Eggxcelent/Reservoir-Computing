#ifndef UTILITY_H
#define UTILITY_H
#define _USE_MATH_DEFINES
#include<iostream>
#include <cmath>
#include<vector>
#include<algorithm>
#include <random>

#include "Eigen/QR"
#include "Eigen/SVD"
#include "Eigen/LU"

namespace  Utility
{
    extern std::random_device rand_dev;
    extern std::mt19937 generator;

    //Sorting three input numbers.
    void Sort(int& a, int& b, int& c);

    //Sorting two input numbers
    void Sort(int& a, int& b);

    //Euclidean distance between two points on x y plane, will overload for 3D
    double Eucl_Dist(double x1, double y1, double x2, double y2);

    //Angle for line between two points.
    double Angle(double x0, double x1, double y0, double y1);

    //Randomly chosen number between N and M
    double Uniform(double M, double N);

    //Randomy chosen log to the base 10 uniform number between initial and final values.
    double Log_10_Uniform(double initial, double finalvalue);


    //X component of vector from one point to another
    double X_Comp(double vectorsum, double theta);

    //Y component of vector from one point to another.
    double Y_Comp(double vectorsum, double theta);


    double Rand_In_Range_Exp(double& var_min, double& var_max);

    //Mean Squared Error between vector A and estimator Ahat
    double MSE(std::vector<double>& A, std::vector<double>& Ahat);

    //Remove duplicates from two dimensional vector
    void Remove_Duplicates(std::vector<std::pair<double,double>>& x); 

    double White_Noise_Generator(bool feedback);
}


#endif 
