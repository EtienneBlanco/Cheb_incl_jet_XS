/*=====================================================\
|                                                      |
|     Functions related to Chebyshev decomposition     |
|	                                                   |
\=====================================================*/


// On this file are gathered all the definitions of the mathematical functions used.

#ifndef CHEBYSHEV_H
#define CHEBYSHEV_H

#include <math.h>

//////////////////////////////////////////////////
//prototype of the used functions               //
double cheb_T(int n, double x);                 //
double cheb_int(int i);                         //
double y_(double x, double a, double b);        //
double y_inv(double x, double a, double b);     //
double y_log(double x, double a, double b);     //
double y_log_inv(double x, double a, double b); //
double cheb_node(int N, int i);                 //
//////////////////////////////////////////////////



//return the value of the Chebyshev polynomial T_n on x
double cheb_T(int n, double x)
{
    double T;
    if (x<-1)
    {
        cout << "\nError cheb_T cannot be evaluated on " << x << endl;
        return 0;
    }
    else if (x>1)
    {
        cout << "\nError cheb_T cannot be evaluated on " << x << endl;
        return 0;
    }
    T = cos(n*acos(x));

    return T;
}

//Return the integral of T_i (over [-1,1])
double cheb_int(int i)
{
    double y;
    if (i==1)
        //y = .5;
        y = 0;
    else
        //y = (i*sin(i*M_PI/2)-1)/double(i*i-1);
        y = (pow(-1., double(i))+1)/double(1-i*i);
    return y;
}

//roots of the Chebyshev polynomials
double cheb_node(int N, int i)
{
    double y = cos(M_PI/N*(i+.5));
    return y;
}

//linear bijection going from [a,b] to [-1,1]
double y_(double x, double a, double b)
{
    double y=(2*x-b-a)/(b-a);
    return y;
}
//inverse of y_
double y_inv(double x, double a, double b)
{
    double y=0.5*((b-a)*x+b+a);
    return y;
}

//logarithmic bijection going from [a,b] to [-1,1]
double y_log(double x, double a, double b)
{
    double y=1+2*log(x/b)/log(b/a);
    return y;
}
//inverse of y_log
double y_log_inv(double x, double a, double b)
{
    double y=b*pow(b/a, .5*(x-1));
    return y;
}

#endif // CHEBYSHEV_H
