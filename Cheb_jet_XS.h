/*====================================================\
|                                                     |
|     Calculation of the cross-section considered     |
|	                                                  |
\====================================================*/


// On this file are gathered the functions related to the calculation of the cross-section.

#ifndef CHEB_XS_H
#define CHEB_XS_H

#include <math.h>
#include "Chebyshev.h"

void cheb_XSy(double *&sigmaAAy, double *&sigmappy, double *&PmeanAA, double *&Pmeanpp, double *&Y, double *&X, double **&factChebX, MyFF FFg2g, MyFF FFg2S, MyFF FFS2S, MyFF FFS2g, LHAPDF::PDF *&PDF1, MyTMD2 TMDg_a, MyTMD2 TMDg_f, parameters param);
void cheb_XSp(double *&sigmaAAp, double *&sigmappp, double *&PT, double *&X, double **&factChebX, MyFF FFg2g, MyFF FFg2S, MyFF FFS2S, MyFF FFS2g, LHAPDF::PDF *&PDF1, MyTMD2 TMDg_a, MyTMD2 TMDg_f, parameters param);

void trap_XSy(double **&sigma, double *&sigmay, double *&PT, parameters param);
void trap_XSp(double **&sigma, double *&sigmap, double *&Y, parameters param);
void stepright_XSy(double **&sigma, double *&sigmay, double *&PT, parameters param);
void stepleft_XSy(double **&sigma, double *&sigmay, double *&PT, parameters param);
void stepright_XSp(double **&sigma, double *&sigmap, double *&Y, parameters param);
void stepleft_XSp(double **&sigma, double *&sigmap, double *&Y, parameters param);

void trap_Pmean(double **&sigma, double *&sigmay, double *&Pmean, double *&PT, parameters param);
void stepright_Pmean(double **&sigma, double *&sigmay, double *&Pmean, double *&PT, parameters param);
void stepleft_Pmean(double **&sigma, double *&sigmay, double *&Pmean, double *&PT, parameters param);

double tau(double qt, double y, parameters param); // Proper time of the Fragmentation Functions (alphabar from parameters)
double partonicXS_gjet(double qt, double yq, LHAPDF::PDF *&PDFg, MyTMD2 TMDg_a, parameters param); // Partonic cross-section for gluon jet
double partonicXS_qjet(double qt, double yq, LHAPDF::PDF *&PDF1, MyTMD2 TMDg_f, int Nf, parameters param); // Partonic cross-section for quark jet
double M2gsg_g(double qt2, parameters param); //Matrix elements squared for g*g->g
double M2gsq_q(double qt2, parameters param); //Matrix elements squared for g*q->q

void cheb_XS(double **&sigmaAA, double *&sigmaAAy, double *&sigmaAAp, double **&sigmapp, double *&sigmappy, double *&sigmappp, double *&PmeanAA, double *&Pmeanpp, double *&Y, double *&PT, MyFF FFg2g, MyFF FFg2S, MyFF FFS2S, MyFF FFS2g, LHAPDF::PDF *&PDF1, MyTMD2 TMDg_a, MyTMD2 TMDg_f, parameters param)
{
    const int Nx = param.Nx;
    const int Ny = param.Ny;
    const int Np = param.Np;
    const int Nf = param.Nf;

    double eps = param.eps;
    double qmin = param.qtmin;
    double qmax = param.qtmax;

    //Memory allocation
    double *X = new double[Nx]; //Values in which x is evaluated (inverse image of the Chebyshev nodes)
    double T[Nx][Nx]; //Chebyshev polynomials T_i,j=T_i(y_j)
    double **factChebX = new double*[Nx]; //Coefficient depending only on i,j

    double ***sigma_part_g = new double**[Ny]; //used points of the partonic cross section for gluon jet
    double ***sigma_part_S = new double**[Ny]; //used points of the partonic cross section for quark Singlet jet

    double ***D_g = new double**[Ny]; //Point of the gluon fragmentation function needed
    double ***D_S = new double**[Ny]; //Point of the quark Singlet fragmentation function needed

    //Calculation of recurrent quantities
    for (int i = 0; i < Nx; ++i)
    {
        factChebX[i] = new double[Nx];

        if (eps == 0)
            X[i] = y_inv(cheb_node(Nx, i),0,1);
        else
            X[i] = y_log_inv(cheb_node(Nx, i),eps,1);

        for (int j = 0; j < Nx; ++j)
        {
            // /!\ /!\ /!\
            //Can be improved by iterative definition of Chebyshev polynomials
            T[i][j] = cheb_T(i, cheb_node(Nx,j));

            if (eps == 0)
                factChebX[i][j] = T[i][j]*cheb_int(i)/Nx;
            else
                factChebX[i][j] = -T[i][j]*pow(eps,(1-cheb_node(Nx, j))/2.)*log(eps)*cheb_int(i)/Nx;
        }
    }

    //Calculation of dsigmay/dy/dpt
    for (int iy = 0; iy < Ny; ++iy)
    {
        for (int ip = 0; ip < Np; ++ip)
        {
            sigmapp[iy][ip] = partonicXS_gjet(PT[ip], Y[iy], PDF1, TMDg_a, param);
            //sigmapp[iy][ip] += partonicXS_qjet(PT[ip], Y[iy], PDF1, TMDg_f, Nf, param);
        }
    }
    cout << "-dsigma_pp/dy/dpt calculated" << endl;

    //Calculation of dsigmaAA/dy/dpt
    for (int iy = 0; iy < Ny; ++iy)
    {
        //cout << "### ### ### iy = " << iy << endl;

        sigma_part_g[iy] = new double*[Np];
        sigma_part_S[iy] = new double*[Np];

        D_g[iy] = new double*[Np];
        D_S[iy] = new double*[Np];

        for (int ip = 0; ip < Np; ++ip)
        {
            //cout << "### ### ip = " << ip << endl;
            D_g[iy][ip] = new double[Nx];
            D_S[iy][ip] = new double[Nx];

            sigma_part_g[iy][ip] = new double[Nx];
            sigma_part_S[iy][ip] = new double[Nx];
            sigmaAA[iy][ip] = 0;
            //sigmapp[iy][ip] = 0;

            for (int j = 0; j < Nx; ++j)
            {
                //cout << "### j = " << j << endl;
                if (((PT[ip]/X[j] >= qmin) && (PT[ip]/X[j] <= qmax)) || ((PT[ip]/X[j] >= qmin) && (qmax == 0)))
                {
                    sigma_part_g[iy][ip][j] = partonicXS_gjet(PT[ip]/X[j], Y[iy], PDF1, TMDg_a, param);
                    sigma_part_S[iy][ip][j] = partonicXS_qjet(PT[ip]/X[j], Y[iy], PDF1, TMDg_f, Nf, param);
                }
                else
                {
                    sigma_part_g[iy][ip][j] = 0;
                    sigma_part_S[iy][ip][j] = 0;
                }

                D_g[iy][ip][j] = FFg2g.interp(X[j],tau(PT[ip]/X[j], Y[iy], param), param);
                D_g[iy][ip][j] += FFg2S.interp(X[j],tau(PT[ip]/X[j], Y[iy], param), param)*2*param.Nf;
                //cout << "D_g(x=" << X[j] << ",tau=" << tau(PT[ip]/X[j], Y[iy], param) << ") = " << D_g[iy][ip][j] << endl;

                D_S[iy][ip][j] = FFS2S.interp(X[j],tau(PT[ip]/X[j], Y[iy], param), param)*2*param.Nf;
                D_S[iy][ip][j] += FFS2g.interp(X[j],tau(PT[ip]/X[j], Y[iy], param), param)*2*param.Nf;
                //cout << "D_S(x=" << X[j] << ",tau=" << tau(PT[ip]/X[j], Y[iy], param) << ") = " << D_S[iy][ip][j] << endl;

                for (int i = 0; i < Nx; ++i)
                {
                    //cout << "i = " << i << endl;
                    sigmaAA[iy][ip] += 1./X[j]/X[j]*sigma_part_g[iy][ip][j]*D_g[iy][ip][j]*factChebX[i][j];
                    sigmaAA[iy][ip] += 1./X[j]/X[j]*sigma_part_S[iy][ip][j]*D_S[iy][ip][j]*factChebX[i][j];
                }
                sigmaAA[iy][ip] -= .5/X[j]/X[j]*sigma_part_g[iy][ip][j]*D_g[iy][ip][j]*factChebX[0][j];
                sigmaAA[iy][ip] -= .5/X[j]/X[j]*sigma_part_S[iy][ip][j]*D_S[iy][ip][j]*factChebX[0][j];
            }
        }
    }
    cout << "-dsigma_AA/dy/dpt calculated" << endl;

    int intmethod = param.intmethod;
    //Integration through trapezoidal method
    if (intmethod == 0)
    {
        trap_XSy(sigmaAA, sigmaAAy, PT, param);
        cout << "-dsigma_AA/dy calculated" << endl;
        trap_XSy(sigmapp, sigmappy, PT, param);
        cout << "-dsigma_pp/dy calculated" << endl;
        trap_XSp(sigmaAA, sigmaAAp, Y, param);
        cout << "-dsigma_AA/dpt calculated" << endl;
        trap_XSp(sigmapp, sigmappp, Y, param);
        cout << "-dsigma_pp/dpt calculated" << endl;

        trap_Pmean(sigmaAA, sigmaAAy, PmeanAA, PT, param);
        cout << "-<pt>_AA calculated" << endl;
        trap_Pmean(sigmapp, sigmappy, Pmeanpp, PT, param);
        cout << "-<pt>_pp calculated" << endl;
    }
    //Integration through Gauss-Chebyshev method
    else if (intmethod == 1)
    {
        cheb_XSy(sigmaAAy, sigmappy, PmeanAA, Pmeanpp, Y, X, factChebX, FFg2g, FFg2S, FFS2S, FFS2g, PDF1, TMDg_a, TMDg_f, param);
        cheb_XSp(sigmaAAp, sigmappp, PT, X, factChebX, FFg2g, FFg2S, FFS2S, FFS2g, PDF1, TMDg_a, TMDg_f, param);
    }
    //Left Riemann sum integration
    else if (intmethod == 2)
    {
        stepleft_XSy(sigmaAA, sigmaAAy, PT, param);
        cout << "-dsigma_AA/dy calculated" << endl;
        stepleft_XSy(sigmapp, sigmappy, PT, param);
        cout << "-dsigma_pp/dy calculated" << endl;
        stepleft_XSy(sigmaAA, sigmaAAp, Y, param);
        cout << "-dsigma_AA/dpt calculated" << endl;
        stepleft_XSy(sigmapp, sigmappp, Y, param);
        cout << "-dsigma_pp/dpt calculated" << endl;

        stepleft_Pmean(sigmaAA, sigmaAAy, PmeanAA, PT, param);
        cout << "-<pt>_AA calculated" << endl;
        stepleft_Pmean(sigmapp, sigmappy, Pmeanpp, PT, param);
        cout << "-<pt>_pp calculated" << endl;
    }
    //Right Riemann sum integration
    else if (intmethod == 3)
    {
        stepright_XSy(sigmaAA, sigmaAAy, PT, param);
        cout << "-dsigma_AA/dy calculated" << endl;
        stepright_XSy(sigmapp, sigmappy, PT, param);
        cout << "-dsigma_pp/dy calculated" << endl;
        stepright_XSy(sigmaAA, sigmaAAp, Y, param);
        cout << "-dsigma_AA/dpt calculated" << endl;
        stepright_XSy(sigmapp, sigmappp, Y, param);
        cout << "-dsigma_pp/dpt calculated" << endl;

        stepright_Pmean(sigmaAA, sigmaAAy, PmeanAA, PT, param);
        cout << "-<pt>_AA calculated" << endl;
        stepright_Pmean(sigmapp, sigmappy, Pmeanpp, PT, param);
        cout << "-<pt>_pp calculated" << endl;
    }

    //*
    for (int iy = 0; iy < Ny; ++iy)
    {
        for (int ip = 0; ip < Np; ++ip)
        {
            delete[] sigma_part_g[iy][ip];
            delete[] sigma_part_S[iy][ip];
            delete[] D_g[iy][ip];
            delete[] D_S[iy][ip];
        }
        delete[] sigma_part_g[iy];
        delete[] sigma_part_S[iy];
        delete[] D_g[iy];
        delete[] D_S[iy];
    }
    delete[] sigma_part_g;
    delete[] sigma_part_S;
    delete[] D_g;
    delete[] D_S;//*/

    return;
}

//Calculation of the differential cross-section dsigma/dy. The integration is done using trapezoid method.
void trap_XSy(double **&sigma, double *&sigmay, double *&PT, parameters param)
{
    const int Ny = param.Ny;
    const int Np = param.Np;

    for (int i = 0; i < Ny; ++i)
    {
        sigmay[i] = 0;
        for (int j = 0; j < Np-1; ++j)
            sigmay[i] += (PT[j+1]-PT[j])*(sigma[i][j]+sigma[i][j+1])/2.;
    }
}

//Calculation of the differential cross-section dsigma/dy. The integration is done using Right Riemann sum.
void stepright_XSy(double **&sigma, double *&sigmay, double *&PT, parameters param)
{
    const int Ny = param.Ny;
    const int Np = param.Np;

    for (int i = 0; i < Ny; ++i)
    {
        sigmay[i] = 0;
        for (int j = 0; j < Np-1; ++j)
            sigmay[i] += (PT[j+1]-PT[j])*sigma[i][j+1];
    }
}

//Calculation of the differential cross-section dsigma/dy. The integration is done using Left Riemann sum.
void stepleft_XSy(double **&sigma, double *&sigmay, double *&PT, parameters param)
{
    const int Ny = param.Ny;
    const int Np = param.Np;

    for (int i = 0; i < Ny; ++i)
    {
        sigmay[i] = 0;
        for (int j = 0; j < Np-1; ++j)
            sigmay[i] += (PT[j+1]-PT[j])*sigma[i][j];
    }
}

//Calculation of the differential cross-section dsigma/dpt. The integration is done using trapezoid method.
void trap_XSp(double **&sigma, double *&sigmap, double *&Y, parameters param)
{
    const int Ny = param.Ny;
    const int Np = param.Np;

    for (int j = 0; j < Np; ++j)
    {
        sigmap[j] = 0;
        for (int i = 0; i < Ny-1; ++i)
            sigmap[j] += (Y[i+1]-Y[i])*(sigma[i][j]+sigma[i+1][j])/2.;
    }
}

//Calculation of the differential cross-section dsigma/dpt. The integration is done using Left Riemann sum.
void stepleft_XSp(double **&sigma, double *&sigmap, double *&Y, parameters param)
{
    const int Ny = param.Ny;
    const int Np = param.Np;

    for (int j = 0; j < Np; ++j)
    {
        sigmap[j] = 0;
        for (int i = 0; i < Ny-1; ++i)
            sigmap[j] += (Y[i+1]-Y[i])*sigma[i][j];
    }
}

//Calculation of the differential cross-section dsigma/dpt. The integration is done using Right Riemann sum.
void stepright_XSp(double **&sigma, double *&sigmap, double *&Y, parameters param)
{
    const int Ny = param.Ny;
    const int Np = param.Np;

    for (int j = 0; j < Np; ++j)
    {
        sigmap[j] = 0;
        for (int i = 0; i < Ny-1; ++i)
            sigmap[j] += (Y[i+1]-Y[i])*sigma[i+1][j];
    }
}

//Calculation of the differential cross-section dsigma/dy. The integration is done using trapezoid method.
void trap_Pmean(double **&sigma, double *&sigmay, double *&Pmean, double *&PT, parameters param)
{
    const int Ny = param.Ny;
    const int Np = param.Np;

    for (int i = 0; i < Ny; ++i)
    {
        Pmean[i] = 0;
        for (int j = 0; j < Np-1; ++j)
            Pmean[i] += (PT[j+1]+PT[j])/2.*(PT[j+1]*sigma[i][j]+PT[j]*sigma[i][j+1])/2./sigmay[i];
    }
}

//Calculation of the differential cross-section dsigma/dy. The integration is done using Right Riemann sum.
void stepright_Pmean(double **&sigma, double *&sigmay, double *&Pmean, double *&PT, parameters param)
{
    const int Ny = param.Ny;
    const int Np = param.Np;

    for (int i = 0; i < Ny; ++i)
    {

        Pmean[i] = 0;
        for (int j = 0; j < Np-1; ++j)
            Pmean[i] += (PT[j+1]+PT[j])/2.*PT[j+1]*sigma[i][j+1]/sigmay[i];
    }
}

//Calculation of the differential cross-section dsigma/dy. The integration is done using Left Riemann sum.
void stepleft_Pmean(double **&sigma, double *&sigmay, double *&Pmean, double *&PT, parameters param)
{
    const int Ny = param.Ny;
    const int Np = param.Np;

    for (int i = 0; i < Ny; ++i)
    {
        Pmean[i] = 0;
        for (int j = 0; j < Np-1; ++j)
            Pmean[i] += (PT[j+1]+PT[j])/2.*PT[j]*sigma[i][j]/sigmay[i];
    }
}

//Calculation of the differential cross-section dsigma/dy and <p>. The integration is done using Chebyshev polynomials expansion.
void cheb_XSy(double *&sigmaAAy, double *&sigmappy, double *&PmeanAA, double *&Pmeanpp, double *&Y, double *&X, double **&factChebX, MyFF FFg2g, MyFF FFg2S, MyFF FFS2S, MyFF FFS2g, LHAPDF::PDF *&PDF1, MyTMD2 TMDg_a, MyTMD2 TMDg_f, parameters param)
{
    const int Nx = param.Nx;
    const int Ny = param.Ny;
    const int Np = param.Np;
    const int Nf = param.Nf;

    double pmin = param.ptmin;
    double pmax = param.ptmax;

    //Memory allocation
    double PT[Np]; //Values in which pt is evaluated for the integral (inverse image of the Chebyshev nodes)
    double T[Np][Np]; //Chebyshev polynomials T_i,j=T_i(y_j)
    double factChebPT[Np][Np]; //Coefficient depending only on i,j

    double ***sigma_part_g = new double**[Ny]; //used points of the partonic cross section for gluon jet
    double ***sigma_part_S = new double**[Ny]; //used points of the partonic cross section for quark Singlet jet

    double ***D_g = new double**[Ny]; //Point of the gluon fragmentation function needed
    double ***D_S = new double**[Ny]; //Point of the quark Singlet fragmentation function needed

    double v; //factor in the Chebyshev expansion

    //Calculation of recurrent quantities
    for (int i = 0; i < Np; ++i)
    {
        PT[i] = y_inv(cheb_node(Np, i), pmin, pmax);

        for (int j = 0; j < Np; ++j)
        {
            T[i][j] = cheb_T(i, cheb_node(Np,j));
            factChebPT[i][j] = (pmax-pmin)*T[i][j]*cheb_int(i)/Np;
        }
    }

    //Calculation of dsigmapp/dy
    for (int iy = 0; iy < Ny; ++iy)
    {
        sigmappy[iy] = 0;
        Pmeanpp[iy] = 0;
        for (int i = 0; i < Np; ++i)
        {
            for (int j = 0; j < Np; ++j)
            {
                sigmappy[iy] += partonicXS_gjet(PT[j], Y[iy], PDF1, TMDg_a, param)*factChebPT[i][j];
                sigmappy[iy] += partonicXS_qjet(PT[j], Y[iy], PDF1, TMDg_f, Nf, param)*factChebPT[i][j];

                Pmeanpp[iy] += PT[i]*partonicXS_gjet(PT[j], Y[iy], PDF1, TMDg_a, param)*factChebPT[i][j];
                Pmeanpp[iy] += PT[i]*partonicXS_qjet(PT[j], Y[iy], PDF1, TMDg_f, Nf, param)*factChebPT[i][j];
            }
            sigmappy[iy] -= .5*partonicXS_gjet(PT[i], Y[iy], PDF1, TMDg_a, param)*factChebPT[0][i];
            sigmappy[iy] -= .5*partonicXS_qjet(PT[i], Y[iy], PDF1, TMDg_f, Nf, param)*factChebPT[0][i];

            Pmeanpp[iy] -= .5*PT[i]*partonicXS_gjet(PT[i], Y[iy], PDF1, TMDg_a, param)*factChebPT[0][i];
            Pmeanpp[iy] -= .5*PT[i]*partonicXS_qjet(PT[i], Y[iy], PDF1, TMDg_f, Nf, param)*factChebPT[0][i];
        }
    }
    cout << "-dsigma_pp/dy and <pt>_pp calculated" << endl;

    //Calculation of dsigmaAA/dy
    for (int iy = 0; iy < Ny; ++iy)
    {
        sigma_part_g[iy] = new double*[Np];
        sigma_part_S[iy] = new double*[Np];
        D_g[iy] = new double*[Np];
        D_S[iy] = new double*[Np];
        sigmaAAy[iy] = 0;
        PmeanAA[iy] = 0;

        for (int jp = 0; jp < Np; ++jp)
        {
            sigma_part_g[iy][jp] = new double[Nx];
            sigma_part_S[iy][jp] = new double[Nx];
            D_g[iy][jp] = new double[Nx];
            D_S[iy][jp] = new double[Nx];

            for (int jx = 0; jx < Nx; ++jx)
            {
                sigma_part_g[iy][jp][jx] = partonicXS_gjet(PT[jp]/X[jx], Y[iy], PDF1, TMDg_a, param);
                sigma_part_S[iy][jp][jx] = partonicXS_qjet(PT[jp]/X[jx], Y[iy], PDF1, TMDg_f, Nf, param);

                D_g[iy][jp][jx] = FFg2g.interp(X[jx],tau(PT[jp]/X[jx], Y[iy], param), param);
                D_g[iy][jp][jx] += FFg2S.interp(X[jx],tau(PT[jp]/X[jx], Y[iy], param), param);

                D_S[iy][jp][jx] = FFS2S.interp(X[jx],tau(PT[jp]/X[jx], Y[iy], param), param);
                D_S[iy][jp][jx] += FFS2g.interp(X[jx],tau(PT[jp]/X[jx], Y[iy], param), param);

                for (int ix = 0; ix < Nx; ++ix)
                {
                    for (int ip = 0; ip < Np; ++ip)
                    {
                        if ((ix == 0) && (ip == 0))
                            v = .25;
                        else if ((ix == 0) || (ip == 0))
                            v = .5;
                        else
                            v = 1.;
                        sigmaAAy[iy] += v/X[jx]/X[jx]*sigma_part_g[iy][jp][jx]*D_g[iy][jp][jx]*factChebX[ix][jx]*factChebPT[ip][jp];
                        sigmaAAy[iy] += v/X[jx]/X[jx]*sigma_part_S[iy][jp][jx]*D_S[iy][jp][jx]*factChebX[ix][jx]*factChebPT[ip][jp];

                        PmeanAA[iy] += v*PT[jp]/X[jx]/X[jx]*sigma_part_g[iy][jp][jx]*D_g[iy][jp][jx]*factChebX[ix][jx]*factChebPT[ip][jp];
                        PmeanAA[iy] += v*PT[jp]/X[jx]/X[jx]*sigma_part_S[iy][jp][jx]*D_S[iy][jp][jx]*factChebX[ix][jx]*factChebPT[ip][jp];
                    }
                }
            }
        }
    }
    cout << "-dsigma_AA/dy and <pt>_AA calculated" << endl;

    for (int iy = 0; iy < Ny; ++iy)
    {
        for (int ip = 0; ip < Np; ++ip)
        {
            delete[] sigma_part_g[iy][ip];
            delete[] sigma_part_S[iy][ip];
            delete[] D_g[iy][ip];
            delete[] D_S[iy][ip];
        }
        delete[] sigma_part_g[iy];
        delete[] sigma_part_S[iy];
        delete[] D_g[iy];
        delete[] D_S[iy];
    }
    delete[] sigma_part_g;
    delete[] sigma_part_S;
    delete[] D_g;
    delete[] D_S;

    return;
}

//Calculation of the differential cross-section dsigma/dpt. The integration is done using Chebyshev polynomials expansion.
void cheb_XSp(double *&sigmaAAp, double *&sigmappp, double *&PT, double *&X, double **&factChebX, MyFF FFg2g, MyFF FFg2S, MyFF FFS2S, MyFF FFS2g, LHAPDF::PDF *&PDF1, MyTMD2 TMDg_a, MyTMD2 TMDg_f, parameters param)
{
    const int Nx = param.Nx;
    const int Ny = param.Ny;
    const int Np = param.Np;
    const int Nf = param.Nf;

    double ymin = param.ymin;
    double ymax = param.ymax;

    //Memory allocation
    double Y[Ny]; //Values in which pt is evaluated for the integral (inverse image of the Chebyshev nodes)
    double T[Ny][Ny]; //Chebyshev polynomials T_i,j=T_i(y_j)
    double factChebY[Ny][Ny]; //Coefficient depending only on i,j

    double ***sigma_part_g = new double**[Np]; //used points of the partonic cross section for gluon jet
    double ***sigma_part_S = new double**[Np]; //used points of the partonic cross section for quark Singlet jet

    double ***D_g = new double**[Np]; //Point of the gluon fragmentation function needed
    double ***D_S = new double**[Np]; //Point of the quark Singlet fragmentation function needed

    double v; //factor in the Chebyshev expansion

    //Calculation of recurrent quantities
    for (int i = 0; i < Ny; ++i)
    {
        Y[i] = y_inv(cheb_node(Ny, i), ymin, ymax);

        for (int j = 0; j < Ny; ++j)
        {
            T[i][j] = cheb_T(i, cheb_node(Ny,j));
            factChebY[i][j] = (ymax-ymin)*T[i][j]*cheb_int(i)/Ny;
        }
    }

    //Calculation of dsigmapp/dp
    for (int ip = 0; ip < Np; ++ip)
    {
        sigmappp[ip] = 0;
        for (int i = 0; i < Ny; ++i)
        {
            for (int j = 0; j < Ny; ++j)
            {
                sigmappp[ip] += partonicXS_gjet(PT[ip], Y[j], PDF1, TMDg_a, param)*factChebY[i][j];
                sigmappp[ip] += partonicXS_qjet(PT[ip], Y[j], PDF1, TMDg_f, Nf, param)*factChebY[i][j];
            }
            sigmappp[ip] -= .5*partonicXS_gjet(PT[ip], Y[i], PDF1, TMDg_a, param)*factChebY[0][i];
            sigmappp[ip] -= .5*partonicXS_qjet(PT[ip], Y[i], PDF1, TMDg_f, Nf, param)*factChebY[0][i];
        }
    }
    cout << "-dsigma_pp/dpt calculated" << endl;

    //Calculation of dsigmaAA/dp
    for (int ip = 0; ip < Np; ++ip)
    {
        sigma_part_g[ip] = new double*[Ny];
        sigma_part_S[ip] = new double*[Ny];
        D_g[ip] = new double*[Ny];
        D_S[ip] = new double*[Ny];
        sigmaAAp[ip] = 0;

        for (int jy = 0; jy < Ny; ++jy)
        {
            sigma_part_g[ip][jy] = new double[Nx];
            sigma_part_S[ip][jy] = new double[Nx];
            D_g[ip][jy] = new double[Nx];
            D_S[ip][jy] = new double[Nx];

            for (int jx = 0; jx < Nx; ++jx)
            {
                sigma_part_g[ip][jy][jx] = partonicXS_gjet(PT[ip]/X[jx], Y[jy], PDF1, TMDg_a, param);
                sigma_part_S[ip][jy][jx] = partonicXS_qjet(PT[ip]/X[jx], Y[jy], PDF1, TMDg_f, Nf, param);

                D_g[ip][jy][jx] = FFg2g.interp(X[jx],tau(PT[ip]/X[jx], Y[jy], param), param);
                D_g[ip][jy][jx] += FFg2S.interp(X[jx],tau(PT[ip]/X[jx], Y[jy], param), param);

                D_S[ip][jy][jx] = FFS2S.interp(X[jx],tau(PT[ip]/X[jx], Y[jy], param), param);
                D_S[ip][jy][jx] += FFS2g.interp(X[jx],tau(PT[ip]/X[jx], Y[jy], param), param);

                for (int ix = 0; ix < Nx; ++ix)
                {
                    for (int iy = 0; iy < Ny; ++iy)
                    {
                        if ((ix == 0) && (iy == 0))
                            v = .25;
                        else if ((ix == 0) || (iy == 0))
                            v = .5;
                        else
                            v = 1.;
                        sigmaAAp[ip] += v/X[jx]/X[jx]*sigma_part_g[ip][jy][jx]*D_g[ip][jy][jx]*factChebX[ix][jx]*factChebY[iy][jy];
                        sigmaAAp[ip] += v/X[jx]/X[jx]*sigma_part_S[ip][jy][jx]*D_S[ip][jy][jx]*factChebX[ix][jx]*factChebY[iy][jy];
                    }
                }
            }
        }
    }
    cout << "-dsigma_AA/dpt calculated" << endl;

    for (int ip = 0; ip < Np; ++ip)
    {
        for (int iy = 0; iy < Ny; ++iy)
        {
            delete[] sigma_part_g[ip][iy];
            delete[] sigma_part_S[ip][iy];
            delete[] D_g[ip][iy];
            delete[] D_S[ip][iy];
        }
        delete[] sigma_part_g[ip];
        delete[] sigma_part_S[ip];
        delete[] D_g[ip];
        delete[] D_S[ip];
    }
    delete[] sigma_part_g;
    delete[] sigma_part_S;
    delete[] D_g;
    delete[] D_S;

    return;
}

// Proper time as a function of qt, y and t (fixed in parameters)
double tau(double qt, double yq, parameters param)
{
    double alphab = param.alphabar;
    double t = param.t;
    double q = param.qhat/param.fmtoGeV; //[GeV/fm2]

    return alphab*t*sqrt(q/(qt*cosh(yq)));
}

//Calculation of the partonic cross-section (for gluon jet production)
double partonicXS_gjet(double qt, double yq, LHAPDF::PDF *&PDF1, MyTMD2 TMDg_a, parameters param)
{
    double x1 = qt*exp(yq)/param.S; //Bjorken x for the PDF
    double x2 = qt*exp(-yq)/param.S; //Bjorken x for the TMD

    double qt2 = qt*qt;

    //cout << "\tx1 = " << x1 << ",\tx2 = " << x2 << endl;
    if ((x1 < 0)||(x2 < 0)||(x1 > 1)||(x2 > 1))
        return 0;
    //if ((x1 < 0)||(x2 < 1.19*1e-8)||(x1 > 1)||(x2 > 0.01))


    double x1f = PDF1->xfxQ2(21, x1, qt2); //gluon PDf at the needed value
    //double f = PDF1.interp(x1, qt2, param); //gluon PDf at the needed value
    //cout << "xf(x=" << x1 << ", mu2=" << qt2 << ") = "<< x1f << endl;

    //double F = TMDg_a.interp(x2, qt, qt2); //gluon TMD at the needed value
    double F = TMDg_a.interp(x2, qt2, param); //gluon TMD at the needed value
    //cout << "F(x=" << x2 << ",qt=" << qt << ") = " << F << endl;

    double alphaS = PDF1->alphasQ2(qt2);
    //double alphaS = 0.3;
    //cout << "aslpha_S = " << PDF1->alphasQ2(qt2) << endl;

    return M_PI/2./qt/qt2*alphaS*M2gsg_g(qt2, param)*x1f*F;
}

//Calculation of the partonic cross-section (for quark jet production)
double partonicXS_qjet(double qt, double yq, LHAPDF::PDF *&PDF1, MyTMD2 TMDg_f, int Nf, parameters param)
{
    double x1 = qt*exp(yq)/param.S; //Bjorken x for the PDF
    double x2 = qt*exp(-yq)/param.S; //Bjorken x for the TMD

    double qt2 = qt*qt;

    //cout << "\tx1 = " << x1 << ",\tx2 = " << x2 << endl;
    if ((x1 < 0)||(x2 < 0)||(x1 > 1)||(x2 > 1))
        return 0;

    double x1f = 0;
    for (int f = 1; f <= Nf; ++f)
        x1f += PDF1->xfxQ2(f, x1, qt2) + PDF1->xfxQ2(-f, x1, qt2); //quark singlet PDf of flavor f at the needed value
    //cout << "xf(x=" << x1 << ", mu2=" << qt2 << ") = "<< x1f << endl;

    //double F = TMDg_f.interp(x2, qt, qt2); //gluon TMD at the needed value
    double F = TMDg_f.interp(x2, qt2, param); //gluon TMD at the needed value
    //cout << "F(x=" << x2 << ",qt=" << qt << ") = " << F << endl;

    double alphaS = PDF1->alphasQ2(qt2);
    //cout << "aslpha_S = " << PDFg->alphasQ2(qt2) << endl;

    return M_PI/2./qt/qt2*alphaS*M2gsq_q(qt2, param)*x1f*F;
}

//Matrix elements squared for g*g->g
double M2gsg_g(double qt2, parameters param)
{
    int Nc = param.Nc;
    return 16*M_PI*Nc/(Nc*Nc-1)*qt2;
}

//Matrix elements squared for g*q->q
double M2gsq_q(double qt2, parameters param)
{
    int Nc = param.Nc;
    return 16*M_PI*param.Cf/(Nc*Nc-1)*qt2;
}

#endif // CHEB_XS_H
