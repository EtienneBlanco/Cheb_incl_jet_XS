/*=============================================\
|                                              |
|  Calculation of inclusive jet cross-section  |
|	                                           |
\=============================================*/



//This program calculate the single inclusive jet cross-section within hybrid factorization for nuclear collision.
//The program is convoluting a given TMD, a given PDF with the partonic cross-section and finally with a given fragmentation function
//(solution of the BDIM equations).
//The partonic cross-section is based on the sum of the squared matrix elements for the following processes : g*g->g, g*q->q

//Compile with : g++ main.cpp -o jet `lhapdf-config --cflags --ldflags` `TMDlib-config --cppflags --ldflags`
//assuming that LHAPDF is installed (tested with LHAPDF 6.4)


#include <fstream>
#include <iostream>
#include <string>
#include  <string.h>
using namespace std;

//LHAPDF
#include "LHAPDF/LHAPDF.h"
//using namespace LHAPDF;

//TMDlib
//#include "tmdlib/TMDlib.h"

//Associated files
#include "Cheb_jet_Parameters.h"
#include "Cheb_jet_Distribution_classes.h"
#include "Cheb_jet_XS.h"
#include "Cheb_jet_Result_grid.h"

int program(string parametersfile); //main program
int help(void); //help display

//argument management
int main(int argc, char** argv)
{
    if (argc == 1)
        //launch the program with default value for the name of the parameters file
        program("Parameters");
    else if (argc == 2)
    {
        //if (argv[1] == "help")
        if (strcmp(argv[1], "help")==0)
            //display help
            help();
        else
        {
            //launch the program with specified name of the parameters file
            return program(argv[1]);
        }
    }
    else
    {
        //wrong input display
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        cout << "!!  /!\\ Wrong argument(s)  !!" << endl;
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

        cout << "\nThis program should be launched with either :" << endl;
        cout << "-no argument (parameters are initialized from the file 'Parameters') :" << endl;
        cout << "\t" << argv[0] << endl;
        cout << "-the name of the parameters file to be use in argument :" << endl;
        cout << "\t" << argv[0] << " 'parameterfile'" << endl;
        cout << "-help argument to display help :" << endl;
        cout << "\t" << argv[0] << " help" << endl;
    }
    return 1;
}

//help
int help(void)
{
    cout << "\n#////////////#" << endl;
    cout << "#//  Help  //#" << endl;
    cout << "#/////////////////////////////////////////////#" << endl;

    cout << "# This program calculate the single inclusive jet cross-section within hybrid factorization for nuclear collision." << endl;
    cout << "# It does the convolution between the given TMDs (for gluons in both fundamental and adjoint representation, either grids or managed by TMDlib)," << endl;
    cout << "# PDF (either grids or managed with LHAPDF), fragmentation functions (as grids), and the hard cross section (where only g*g->g and g*q->q are considered)." << endl;
    cout << "#" << endl;
    cout << "# More details on the grids to use and on the parameters are given in the file 'Parameters' and in the manual." << endl;
    cout << "#" << endl;
    cout << "# The results are given in 9 files for the differential cross-sections (for pp and AA collisions) and for the nuclear modification factor, with the following suffixes : " << endl;
    cout << "# -'_xsAA' for dsigmaa_AA/dpt/dy as a 3 column grid (y, pt[GeV], dsigma[pb/GeV])" << endl;
    cout << "# -'_xsAA_pt' for dsigmaa_AA/dpt as a 2 column grid (pt[GeV], dsigma[pb/GeV])" << endl;
    cout << "# -'_xsAA' for dsigmaa_AA/dy as a 2 column grid (y, dsigma[pb])" << endl;
    cout << "# -'_xspp' for dsigmaa_pp/dpt/dy as a 3 column grid (y, pt[GeV], dsigma[pb/GeV])" << endl;
    cout << "# -'_xspp_pt' for dsigmaa_pp/dpt as a 2 column grid (pt[GeV], dsigma[pb/GeV])" << endl;
    cout << "# -'_xspp' for dsigmaa_pp/dy as a 2 column grid (y, dsigma[pb])" << endl;
    cout << "# -'_RAA' for R_AA/dpt/dy as a 3 column grid (y, pt[GeV], R_AA[GeV^-1])" << endl;
    cout << "# -'_RAA_pt' for R_AA/dpt as a 2 column grid (pt[GeV], R_AA[GeV^-1])" << endl;
    cout << "# -'_RAA' for R_AA/dy as a 2 column grid (y, R_AA)" << endl;
    cout << "#" << endl;
    cout << "# The program should be launched either with the name of the parameters file in argument or without argument." << endl;
    cout << "# In this case, parameters are read from the file 'Parameters' if present or set to default values" << endl;
    cout << "# (but this should be avoided since the default values for the grid path might leads to no file...)." << endl;

    cout << "#/////////////////////////////////////////////#\n" << endl;
}

//main program
int program(string parametersfile)
{
    //Parameters initialization
    ///////////////////////////
    //cout << "Initialization" << endl;
    parameters param;
    param.from_file(parametersfile);

    int Ny = param.Ny;
    int Np = param.Np;
    double ym = param.ymin;
    double yM = param.ymax;
    double ptm = param.ptmin;
    double ptM = param.ptmax;

    //Memory allocation
    ///////////////////
    //cout << "-Memory allocation" << endl;
    double *Y = new double[Ny]; //rapidity grid
    double *PT = new double[Np]; //Transverse momentum grid
    double *Pmeanpp = new double[Ny]; //Mean transverse momentum (Hard part)
    double *PmeanAA = new double[Ny]; //Mean transverse momentum (AA)
    double **sigmaAA = new double*[Ny]; //Cross-section (AA)
    double *sigmaAAy = new double[Ny]; //Cross-section (AA)
    double *sigmaAAp = new double[Np]; //Cross-section (AA)
    double **sigmapp = new double*[Ny]; //Hard cross-section (pp)
    double *sigmappy = new double[Ny]; //Hard cross-section (pp)
    double *sigmappp = new double[Np]; //Hard cross-section (pp)
    double **RAA = new double*[Ny]; // Nuclear modification factor s_AA/s_pp
    double *RAAy = new double[Ny]; // Nuclear modification factor s_AA/s_pp integrated over pt
    double *RAAp = new double[Np]; // Nuclear modification factor s_AA/s_pp integrated over y

    for (int i = 0; i < Ny; ++i)
    {
        sigmaAA[i] = new double[Np];
        sigmapp[i] = new double[Np];
        RAA[i] = new double[Np];
        Y[i] = ym + i*(yM-ym)/(Ny-1);
    }
    for (int i = 0; i < Np; ++i)
        PT[i] =  ptm + i*(ptM-ptm)/(Np-1);
    cout << "\n*Memory allocated*" << endl;

    //Initialization (loading the grids)
    ////////////////
    //MyPDF PDFg;
    //PDFg.create("PDF_test"); //gluon PDF
    cout << "\n*Loading PDF through LHAPDF*" << endl;
    LHAPDF::PDF* PDF1 = LHAPDF::mkPDF(param.PDFset, 0); //Create PDF class (LHAPDF)

    MyTMD2 TMDg_a, TMDg_f;
    TMDg_a.create_readN(param.TMDsetga); //Adjoint gluon TMD
    TMDg_f.create_readN(param.TMDsetgf); //Fundamental gluon TMD

    MyFF FFg2g, FFg2S, FFS2S, FFS2g;
    FFg2g.create_readN(param.FFsetg2g); //gluon fragmentation function
    FFg2S.create_readN(param.FFsetg2S); //gluon fragmentation function
    FFS2S.create_readN(param.FFsetS2S); //quark singlet fragmentation function
    FFS2g.create_readN(param.FFsetS2g); //quark singlet fragmentation function
    cout << "\n*Distributions loaded*" << endl;

    //TMD TMDlib_g_a, TMDlib_g_f;
    /*
    string TMDlibname = "PB-NLO-HERAI+II-2018-set1";
    int irep = 0;
    TMDlib::TMDinit(TMDlibname, irep);//*/

    //Integration
    /////////////
    cout << "\nIntegration" << endl;
    cheb_XS(sigmaAA, sigmaAAy, sigmaAAp, sigmapp, sigmappy, sigmappp, PmeanAA, Pmeanpp, Y, PT, FFg2g, FFg2S, FFS2S, FFS2g, PDF1, TMDg_a, TMDg_f, param);
    cout << "\n*Cross-section calculated*" << endl;

    //Nuclear modification factor
    /////////////////////////////
    for (int i = 0; i < Ny; ++i)
    {
        for (int j = 0; j < Np; ++j)
            RAA[i][j] = sigmaAA[i][j]/sigmapp[i][j];
        RAAy[i] = sigmaAAy[i]/sigmappy[i];
    }
    for (int j = 0; j < Np; ++j)
        RAAp[j] = sigmaAAp[j]/sigmappp[j];
    cout << "*Nuclear modification factor calculated*" << endl;

    double sigmatot = 0;
    for (int i = 0; i < Ny-1; ++i)
        sigmatot += (Y[i+1]-Y[i])*(sigmaAAy[i]+sigmaAAy[i+1])/2.;
    cout << "\nTotal cross-section (integration of dsigma/dy) : " << sigmatot/param.pbtoGeV << " [pb]" << endl;
    sigmatot = 0;
    for (int i = 0; i < Np-1; ++i)
        sigmatot += (PT[i+1]-PT[i])*(sigmaAAp[i]+sigmaAAp[i+1])/2.;
    cout << "Total cross-section (integration of dsigma/dpt) : " << sigmatot/param.pbtoGeV << " [pb]" << endl;

    //Grid files writing
    ////////////////////
    string sigma_filename = param.gridname;
    cout << "\n#// Result grids :" << endl;
    cout << "#/////////////////" << endl;
    cout << "#" << endl;

    cout << "# *Writing dsigma/dy/dpt on : " << sigma_filename+"_xsAA" << endl;
    write_sigma(sigmaAA, Y, PT, Ny, Np, sigma_filename+"_xsAA", param);
    cout << "# *Writing dsigma/dpt on : " << sigma_filename+"_xsAA_pt" << endl;
    write_sigma_pt(sigmaAAp, PT, Np, sigma_filename+"_xsAA_pt", param);
    cout << "# *Writing dsigma/dy on : " << sigma_filename+"_xsAA_y" << endl;
    write_sigma_y(sigmaAAy, Y, Ny, sigma_filename+"_xsAA_y", param);

    cout << "#" << endl;

    cout << "# *Writing dsigma_pp/dy/dpt on : " << sigma_filename+"_xspp" << endl;
    write_sigma(sigmapp, Y, PT, Ny, Np, sigma_filename+"_xspp", param);
    cout << "# *Writing dsigma_pp/dpt on : " << sigma_filename+"_xspp_pt" << endl;
    write_sigma_pt(sigmappp, PT, Np, sigma_filename+"_xspp_pt", param);
    cout << "# *Writing dsigma_pp/dy on : " << sigma_filename+"_xspp_y" << endl;
    write_sigma_y(sigmappy, Y, Ny, sigma_filename+"_xspp_y", param);

    cout << "#" << endl;

    cout << "# *Writing R_AA(y,pt) on : " << sigma_filename+"_RAA" << endl;
    write_RAA(RAA, Y, PT, Ny, Np, sigma_filename+"_RAA");
    cout << "# *Writing R_AA(pt) on : " << sigma_filename+"_RAA_pt" << endl;
    write_RAA_pt(RAAp, PT, Np, sigma_filename+"_RAA_pt");
    cout << "# *Writing R_AA(y) on : " << sigma_filename+"_RAA_y" << endl;
    write_RAA_y(RAAy, Y, Ny, sigma_filename+"_RAA_y");

    cout << "#" << endl;

    cout << "# *Writing <p>_pp(y) on : " << sigma_filename+"_Pmpp" << endl;
    write_sigma_y(Pmeanpp, Y, Ny, sigma_filename+"_Pmpp", param);
    cout << "# *Writing <p>_AA(y) on : " << sigma_filename+"_PmAA" << endl;
    write_sigma_y(PmeanAA, Y, Ny, sigma_filename+"_PmAA", param);

    //Memory liberation
    ///////////////////
    for (int i = 0; i < Ny; ++i)
    {
        delete[] sigmaAA[i];
        delete[] sigmapp[i];
    }
    delete[] sigmaAA;
    delete[] sigmaAAy;
    delete[] sigmaAAp;
    delete[] sigmapp;
    delete[] sigmappy;
    delete[] sigmappp;
    delete[] Y;
    delete[] PT;

    cout << "\n*Freed memory*" << endl;

    return 0;
}
