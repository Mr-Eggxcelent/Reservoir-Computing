#define _USE_MATH_DEFINES
#include<random>
#include <cmath>
#include<iostream>

#include "Eigen/QR"
#include "Eigen/SVD"
#include "Eigen/LU"

namespace Utility {

	std::random_device rand_dev;
	std::mt19937 generator(rand_dev());

    //Sorting three input numbers.
    void Sort(int& a, int& b, int& c)
    {
        if (a > b)
        {
            int tmp = a;
            a = b;
            b = tmp;
        }
        if (a > c)
        {
            int tmp = a;
            a = c;
            c = tmp;
        }
        if (b > c)
        {
            int tmp = b;
            b = c;
            c = tmp;
        }

    }
    //Sorting two input numbers
    void Sort(int& a, int& b)
    {
        if (a > b)
        {
            int tmp = a;
            a = b;
            b = tmp;
        }

    }

    //Euclidean distance between two points on x y plane, will overload for 3D
    double Eucl_Dist(double x1, double y1, double x2, double y2)
    {
        return sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
    }

    //Angle for line between two points.
    double Angle(double x0, double x1, double y0, double y1)
    {
        return atan((y1 - y0) / (x1 - x0));
    }

    //Randomly chosen number between N and M
    double Uniform(double M, double N)
    {
        std::uniform_real_distribution<> distr(M, N);
        return distr(generator);
    }


    //Randomy chosen log to the base 10 uniform number between initial and final values.
    double Log_10_Uniform(double initial, double finalvalue)
    {
        return exp(Uniform(initial, finalvalue) / (2.302585093));
    }


    //X component of vector from one point to another
    double X_Comp(double vectorsum, double theta)
    {
        return vectorsum * cos(theta);

    }

    //Y component of vector from one point to another.
    double Y_Comp(double vectorsum, double theta)
    {
        return vectorsum * sin(theta);
    }


    double Rand_In_Range_Exp(double& var_min, double& var_max)
    {
        /*std::cout << "min is:" << var_min << std::endl;
        std::cout << "max is:" << var_max << std::endl;*/

        double log10min = log10(var_min);
        double log10max = log10(var_max);


        double return_value = ((log10max - log10min) * Uniform(0, 1)) + log10min;
        return pow(10, return_value);
    }


    //Mean Squared Error between vector A and estimator Ahat
    double MSE(std::vector<double>& A, std::vector<double>& Ahat)
    {

        double MSEsum = 0;
        double MSE = 0;
        double Total;

        for (unsigned int i = 0; i < A.size(); i++)
        {
            MSEsum += (A[i] - Ahat[i]) * (A[i] - Ahat[i]);
        }

        Total = (double)A.size();

        //cout <<"Inverse total is: " << Total << endl;
        MSE = (1 / Total) * MSEsum;

        return MSE;
    }

    //Remove duplicates from two dimensional vector
    void Remove_Duplicates(std::vector <std::pair<double,double >> & x)
    {
        std::sort(x.begin(), x.end());
        unsigned int i = 1;
        while (i < x.size())
        {
            if (x[i].second == x[i - 1].second && x[i].first == x[i - 1].first)
            {
                x.erase(x.begin() + (i - 1));
                i -= 1;
            }
            i++;
        }

        std::pair<int,int> a;
        i = 0;
        while (i < x.size())
        {
            a = { x[i].first ,x[i].second };

            for (int j = i+1; j < x.size(); j++)
            {
                if (a.first==x[j].second && a.second == x[j].first) {

                    x.erase(x.begin() + (j));
                    break;
                }
            }
            i++;

        }
    }

    double White_Noise_Generator(bool feedback) {

        if (feedback) {
            std::normal_distribution<double> distr(0, 0.000001);
            return distr(generator);
        }
        else {
            return 0;
        }

    }



}