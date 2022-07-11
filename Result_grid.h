/*====================================================\
|                                                     |
|        Functions related to the result grids        |
|	                                                  |
\====================================================*/


// On this file are gathered all the definitions of the mathematical functions used.

#ifndef RESULT_GRID_H
#define RESULT_GRID_H

void write_header(int Ny, int Np, ofstream &grid)
{
    grid << "//Cross-section grid//" << endl;
    grid << "//////////////////////\n" << endl;

    grid << "//This file contain grid describing the differential cross-section for single inclusive jet production in AA collisions.\n" << endl;

    grid << "//Parameters" << endl;
    grid << "//Ny = " << Ny << endl;
    grid << "//Np = " << Np << "\n" << endl;

    grid << "//y, pt [GeV], sigma(y,pt) [pb/GeV]\n" << endl;
}

void write_sigma(double **&sigma, double *&Y, double *&PT, int Ny, int Np, string sigma_filename, parameters param)
{
    ofstream sigma_grid;
    sigma_grid.open(sigma_filename, ios::out);

    //write_header(Ny, Np, sigma_grid);

    for (int i = 0; i < Ny; ++i)
    {
        for (int j = 0; j < Np; ++j)
            sigma_grid << Y[i] << "\t" << PT[j] << "\t" << sigma[i][j]/param.pbtoGeV << endl; //Writing of the solutions on grids
    }
    sigma_grid.close();
}

void write_header_y(int Ny, ofstream &grid)
{
    grid << "//Cross-section grid//" << endl;
    grid << "//////////////////////\n" << endl;

    grid << "//This file contain grid describing the differential cross-section for single inclusive jet production in AA collisions integrated over pt.\n" << endl;

    grid << "//Parameters" << endl;
    grid << "//Ny = " << Ny << "\n" << endl;

    grid << "//y, sigma(y) [pb]\n" << endl;
}

void write_sigma_y(double *&sigmay, double *&Y, int Ny, string sigma_filename, parameters param)
{
    ofstream sigma_grid;
    sigma_grid.open(sigma_filename, ios::out);

    //write_headery(Ny, sigma_grid);

    for (int i = 0; i < Ny; ++i)
            sigma_grid << Y[i] << "\t" << sigmay[i]/param.pbtoGeV << endl; //Writing of the solutions on grids

    sigma_grid.close();
}

void write_header_pt(int Np, ofstream &grid)
{
    grid << "//Cross-section grid//" << endl;
    grid << "//////////////////////\n" << endl;

    grid << "//This file contain grid describing the differential cross-section for single inclusive jet production in AA collisions integrated over y.\n" << endl;

    grid << "//Parameters" << endl;
    grid << "//Np = " << Np << "\n" << endl;

    grid << "//pt [GeV], sigma(pt) [pb/GeV]\n" << endl;
}

void write_sigma_pt(double *&sigmap, double *&PT, int Np, string sigma_filename, parameters param)
{
    ofstream sigma_grid;
    sigma_grid.open(sigma_filename, ios::out);

    //write_header_pt(Np, sigma_grid);

    for (int i = 0; i < Np; ++i)
            sigma_grid << PT[i] << "\t" << sigmap[i]/param.pbtoGeV << endl; //Writing of the solutions on grids

    sigma_grid.close();
}

void write_header_RAA(int Ny, ofstream &grid)
{
    grid << "//Nuclear modification factor grid//" << endl;
    grid << "////////////////////////////////////\n" << endl;

    grid << "//This file contain grid describing the differential nuclear modification factor for single inclusive jet production integrated over pt.\n" << endl;

    grid << "//Parameters" << endl;
    grid << "//Ny = " << Ny << "\n" << endl;

    grid << "//y, RAA(y)\n" << endl;
}

void write_RAA(double **&RAA, double *&Y, double *&PT, int Ny, int Np, string RAA_filename)
{
    ofstream RAA_grid;
    RAA_grid.open(RAA_filename, ios::out);

    //write_header_RAA_pt(Np, RAA_grid);

    for (int i = 0; i < Ny; ++i)
    {
        for (int j = 0; j < Np; ++j)
            RAA_grid << Y[i] << "\t" << PT[j] << "\t" << RAA[i][j] << endl; //Writing of the solutions on grids
    }
    RAA_grid.close();
}

void write_RAA_pt(double *&RAA, double *&PT, int Np, string RAA_filename)
{
    ofstream RAA_grid;
    RAA_grid.open(RAA_filename, ios::out);

    //write_header_RAA_pt(Np, RAA_grid);

    for (int i = 0; i < Np; ++i)
            RAA_grid << PT[i] << "\t" << RAA[i] << endl; //Writing of the solutions on grids

    RAA_grid.close();
}

void write_RAA_y(double *&RAA, double *&Y, int Ny, string RAA_filename)
{
    ofstream RAA_grid;
    RAA_grid.open(RAA_filename, ios::out);

    //write_header_RAA_pt(Np, RAA_grid);

    for (int i = 0; i < Ny; ++i)
            RAA_grid << Y[i] << "\t" << RAA[i] << endl; //Writing of the solutions on grids

    RAA_grid.close();
}

#endif // RESULT_GRID_H
