/*===================================================\
|                                                    |
|                Distribution classes                |
|	                                                 |
\===================================================*/


// On this file are defined the different distribution classes used (for PDF, TMD, FF) and the related functions.

#ifndef DISTRIBUTION_CLASSES_H
#define DISTRIBUTION_CLASSES_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

void findinterpindex(int *i0, int *i1, double x, double *&X); // find i0 and i1 such as X[i0] < x < X[i1], i1 = i0+1 by dichotomy
void findinterpindex_inv(int *i0, int *i1, double x, double *&X); // find i0 and i1 such as X[i1] < x < X[i0], i0 = i1+1 by dichotomy
void findgridN(string filename, int *Nx, int *Ny); // find the number of point x and y in a grid of the form : x y f(x,y)
void findgridN(string filename, int *Nx, int *Ny, int *Nz); // find the number of point x, y and z in a grid of the form : x y z f(x,y,z)



/////////////////////////////////////
// Distribution class definitions  //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Colinear Parton Distribution Function class
class MyPDF
{
    public:
        int Nx, Nmu2; //
        double x_min, x_max, mu2_min, mu2_max;
        double* X;
        double* Mu2;
        double** f;

        void create(string filename); //creation from grid file with header defining Nx and Nmu2
        void create(string filename, int nx, int nmu2); //creation from grid file without header
        void create_readN(string filename); //creation from grid file without header

        double interp(double x, double mu2, parameters param); //interpolation of the PDF f on x, pt, mu2

        //~MyPDF(); //Destructor
};

//Transverse Momentum Dependant PDF class
class MyTMD
{
    public:
        int Nx, Npt, Nmu2; //
        double x_min, x_max, pt_min, pt_max, mu2_min, mu2_max;
        double* X;
        double* Pt;
        double* Mu2;
        double*** F;

        void create(string filename); //creation from grid file with header defining Nx and Nmu2
        void create(string filename, int nx, int npt, int nmu2); //creation from grid file without header
        void create_readN(string filename); //creation from grid file without header

        double interp(double x, double pt, double mu2, parameters param); //interpolation of the TMD F on x, pt, mu2

        //~MyTMD(); //Destructor
};

//Transverse Momentum Dependant PDF (without scale dependence) class
class MyTMD2
{
    public:
        int Nx, Npt; //
        double x_min, x_max, pt_min, pt_max, mu2_min, mu2_max;
        double* X;
        double* Pt;
        double** F;

        void create(string filename); //creation from grid file with header defining Nx
        void create(string filename, int nx, int npt); //creation from grid file without header
        void create_readN(string filename); //creation from grid file without header

        double interp(double x, double pt, parameters param); //interpolation of the TMD F on x, pt

        //~MyTMD2(); //Destructor
};

//Fragmentation Function class
class MyFF
{
    public:
        int Nx, Nt; //
        double x_min, x_max, t_min, t_max, eps;
        double* X;
        double* T;
        double** D;

        void create(string filename); //creation from grid file with header defining Nx and Nmu2
        void create(string filename, int nx, int nt); //creation from grid file without header
        void create_readN(string filename); //creation from grid file without header

        double interp(double x, double t, parameters param); //interpolation of the FF D on x, t
        double interp_cheb(double x, double t, double epsilon); //interpolation of the FF D on x, t using Chebyshev decomposition

        //~MyFF(); //Destructor
};



///////////////////////////////
// Class creation from file  //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Create PDF from file
//The file must include a header
void MyPDF::create(string filename)
{
    cout << "\n\n#//  Loading PDF from file : " << endl;
    cout << "#///////////////////////////" << endl;
    cout << "# Grid : " << filename << endl;

    ifstream pdf_file;
    pdf_file.open(filename);

    string line;
    int i = 0;
    int i_x = 0, i_mu2 = 0;

    if (pdf_file.is_open())
    {
        while (getline(pdf_file, line))
        {
            if (!(line[0]=='/' && line[1]=='/'))
            {
                if (i==0)
                {
                    Nx = strtod(line.c_str(), NULL);
                    i++;
                }
                else if (i==1)
                {
                    Nmu2 = strtod(line.c_str(), NULL);
                    i++;

                    cout << "#\n# Number of points found : " << endl;
                    cout << "# Nx = " << Nx << ", \tNmu2 = " << Nmu2 << ", \tTotal Nx*Nmu2 = " << Nx*Nmu2 << endl;

                    //Memory allocation of the grids in the class PDF
                    f = new double*[Nx];
                    for (int k = 0; k < Nx; ++k)
                        f[k] = new double[Nmu2];

                    X = new double[Nx];
                    Mu2 = new double[Nmu2];
                }
                else
                {
                    istringstream iss(line);
                    iss >> X[i_x] >> Mu2[i_mu2] >> f[i_x][i_mu2];
                    //cout << "\n x = " << X[i_x] << ",\t mu2 = " << Mu2[i_mu2] << ", \t f = " << f[i_x][i_mu2];

                    i_mu2++;
                    if (i_mu2 == Nmu2)
                    {
                        i_mu2=0;
                        i_x++;
                        //if (i_x == Nx)
                        //    i_x=0;
                    }
                }
            }
        }
    }
    pdf_file.close();

    //Variable range :
    x_min = X[0];
    x_max = X[Nx-1];
    mu2_min = Mu2[0];
    mu2_max = Mu2[Nmu2-1];

    cout << "#\n# Variable range : " << endl;
    cout << "# x_min = " << x_min << ",\tx_max = " << x_max << endl;
    cout << "# mu2_min = " << mu2_min << ",\tmu2_max = " << mu2_max << endl;
}

//Create PDF from file
//No header in the file but the number of point in the grid must be given
void MyPDF::create(string filename, int nx, int nmu2)
{
    Nx = nx;
    Nmu2 = nmu2;

    ifstream pdf_file;
    pdf_file.open(filename);

    string line;
    int i_x = 0, i_mu2 = 0;

    //Memory allocation of the grids in the class PDF
    f = new double*[Nx];
    for (int k = 0; k < Nx; ++k)
        f[k] = new double[Nmu2];

    X = new double[Nx];
    Mu2 = new double[Nmu2];

    if (pdf_file.is_open())
    {
        while (getline(pdf_file, line))
        {
            //if (!(line[0]=='/' && line[1]=='/'))
            //{
                istringstream iss(line);
                iss >> X[i_x] >> Mu2[i_mu2] >> f[i_x][i_mu2];
                //cout << "\n x = " << X[i_x] << ",\t mu2 = " << Mu2[i_mu2] << ", \t f = " << f[i_x][i_mu2];

                i_mu2++;
                if (i_mu2 == Nmu2)
                {
                    i_mu2=0;
                    i_x++;
                    //if (i_x == Nx)
                    //    i_x=0;
                }
            //}
        }
    }
    pdf_file.close();

    //Variable range :
    x_min = X[0];
    x_max = X[Nx-1];
    mu2_min = Mu2[0];
    mu2_max = Mu2[Nmu2-1];
}

//Create PDF from file
//Here, the file is read once to find Nx and Nmu2 from its structure
//considering the data (x | mu2 | f(x,mu2)) organized as Nx blocks of size Nmu2 (i.e mu2 increase first at fixed x and reloop when x increase)
//and read a second time to fill X, Mu2 and f
void MyPDF::create_readN(string filename)
{
    cout << "\n\n#//  Loading PDF from file : " << endl;
    cout << "#///////////////////////////" << endl;
    cout << "# Grid : " << filename << endl;

    Nx = 0;
    Nmu2 = 0;

    findgridN(filename, &Nx, &Nmu2); //Determines Nx and Nmu2 from file structure

    cout << "#\n# Number of points found : " << endl;
    cout << "# Nx = " << Nx << ", \tNmu2 = " << Nmu2 << ", \tTotal Nx*Nmu2 = " << Nx*Nmu2 << endl;

    create(filename, Nx, Nmu2); //Copy grid data in the appropriate class variables

    cout << "#\n# Variable range : ";
    cout << "# x_min = " << x_min << ",\tx_max = " << x_max << endl;
    cout << "# mu2_min = " << mu2_min << ",\tmu2_max = " << mu2_max << endl;
}


//Create TMD from file
void MyTMD::create(string filename)
{
    cout << "\n\n#//  Loading TMD from file : " << endl;
    cout << "#///////////////////////////" << endl;
    cout << "# Grid : " << filename << endl;

    ifstream tmd_file;
    tmd_file.open(filename);

    string line;
    int i = 0;
    int i_x = 0, i_pt = 0, i_mu2 = 0;

    if (tmd_file.is_open())
    {
        while (getline(tmd_file, line))
        {
            if (!(line[0]=='/' && line[1]=='/'))
            {
                if (i==0)
                {
                    Nx = strtod(line.c_str(), NULL);
                    i++;
                }
                else if (i==1)
                {
                    Npt = strtod(line.c_str(), NULL);
                    i++;
                }
                else if (i==2)
                {
                    Nmu2 = strtod(line.c_str(), NULL);
                    i++;

                    cout << "#\n# Number of points found : " << endl;
                    cout << "# Nx = " << Nx << ", \tNpt = " << Npt << ", \tNmu2 = " << Nmu2 << ", \tTotal Nx*Npt*Nmu2 = " << Nx*Npt*Nmu2 << endl;

                    //Memory allocation of the grids in the class TMD
                    F = new double**[Nx];
                    for (int j = 0; j < Nx; ++j)
                    {
                        F[j] = new double*[Npt];
                        for (int k = 0; k < Npt; ++k)
                            F[j][k] = new double[Nmu2];
                    }
                    X = new double[Nx];
                    Pt = new double[Npt];
                    Mu2 = new double[Nmu2];
                }
                else
                {
                    istringstream iss(line);
                    iss >> X[i_x] >> Pt[i_pt] >> Mu2[i_mu2] >> F[i_x][i_pt][i_mu2];
                    //cout << "\n x = " << X[i_x] << ",\t pt = " << Pt[i_pt] << ",\t mu2 = " << Mu2[i_mu2] << ", \t F = " << F[i_x][i_pt][i_mu2];

                    i_mu2++;
                    if (i_mu2 == Nmu2)
                    {
                        i_mu2=0;
                        i_pt++;
                        if (i_pt == Npt)
                        {
                            i_pt=0;
                            i_x++;
                            //if (i_x == Nx)
                            //    i_x=0;
                        }
                    }
                }
            }
        }
    }
    tmd_file.close();

    //Variable range :
    x_min = X[0];
    x_max = X[Nx-1];
    pt_min = Pt[0];
    pt_max = Pt[Npt-1];
    mu2_min = Mu2[0];
    mu2_max = Mu2[Nmu2-1];

    cout << "#\n# Variable range : " << endl;
    cout << "# x_min = " << x_min << ",\tx_max = " << x_max << endl;
    cout << "# pt_min = " << pt_min << ",\tpt_max = " << pt_max << endl;
    cout << "# mu2_min = " << mu2_min << ",\tmu2_max = " << mu2_max << endl;
}

//Create TMD from file
//No header in the file but the number of point in the grid must be given
void MyTMD::create(string filename, int nx, int npt, int nmu2)
{
    Nx = nx;
    Npt = npt;
    Nmu2 = nmu2;

    ifstream tmd_file;
    tmd_file.open(filename);

    string line;
    int i_x = 0, i_pt = 0, i_mu2 = 0;

    //Memory allocation of the grids in the class TMD
    F = new double**[Nx];
    for (int j = 0; j < Nx; ++j)
    {
        F[j] = new double*[Npt];
        for (int k = 0; k < Npt; ++k)
            F[j][k] = new double[Nmu2];
    }

    X = new double[Nx];
    Pt = new double[Npt];
    Mu2 = new double[Nmu2];

    if (tmd_file.is_open())
    {
        while (getline(tmd_file, line))
        {
            //if (!(line[0]=='/' && line[1]=='/'))
            //{
                istringstream iss(line);
                iss >> X[i_x] >> Pt[i_pt] >> Mu2[i_mu2] >> F[i_x][i_pt][i_mu2];
                //cout << "\n x = " << X[i_x] << ",\t pt = " << Pt[i_pt] << ",\t mu2 = " << Mu2[i_mu2] << ", \t F = " << F[i_x][i_pt][i_mu2];

                    i_mu2++;
                    if (i_mu2 == Nmu2)
                    {
                        i_mu2=0;
                        i_pt++;
                        if (i_pt == Npt)
                        {
                            i_pt=0;
                            i_x++;
                            //if (i_x == Nx)
                            //    i_x=0;
                        }
                    }
            //}
        }
    }
    tmd_file.close();

    //Variable range :
    x_min = X[0];
    x_max = X[Nx-1];
    pt_min = Pt[0];
    pt_max = Pt[Npt-1];
    mu2_min = Mu2[0];
    mu2_max = Mu2[Nmu2-1];
}

//Create TMD from file
//Here, the file is read once to find Nx, Npt and Nmu2 from its structure (generalize from the PDF routine)
void MyTMD::create_readN(string filename)
{
    cout << "\n\n#//  Loading TMD from file : " << endl;
    cout << "#///////////////////////////" << endl;
    cout << "# Grid : " << filename << endl;

    Nx = 0;
    Npt = 0;
    Nmu2 = 0;

    findgridN(filename, &Nx, &Npt, &Nmu2); //Determines Nx, Npt and Nmu2 from file structure

    cout << "#\n# Number of points found : " << endl;
    cout << "# Nx = " << Nx << ", \tNpt = " << Npt << ", \tNmu2 = " << Nmu2 << ", \tTotal Nx*Npt*Nmu2 = " << Nx*Npt*Nmu2 << endl;

    create(filename, Nx, Npt, Nmu2); //Copy grid data in the appropriate class variables

    cout << "#\n# Variable range : " << endl;
    cout << "# x_min = " << x_min << ",\tx_max = " << x_max << endl;
    cout << "# pt_min = " << pt_min << ",\tpt_max = " << pt_max << endl;
    cout << "# mu2_min = " << mu2_min << ",\tmu2_max = " << mu2_max << endl;
}


//Create TMD from file
void MyTMD2::create(string filename)
{
    cout << "\n\n#//  Loading TMD from file : " << endl;
    cout << "#///////////////////////////" << endl;
    cout << "# Grid : " << filename << endl;

    ifstream tmd_file;
    tmd_file.open(filename);

    string line;
    int i = 0;
    int i_x = 0, i_pt = 0;

    if (tmd_file.is_open())
    {
        while (getline(tmd_file, line))
        {
            if (!(line[0]=='/' && line[1]=='/'))
            {
                if (i==0)
                {
                    Nx = strtod(line.c_str(), NULL);
                    i++;
                }
                else if (i==1)
                {
                    Npt = strtod(line.c_str(), NULL);
                    i++;

                    cout << "#\n# Number of points found : " << endl;
                    cout << "# Nx = " << Nx << ", \tNpt = " << Npt << ", \tTotal Nx*Npt = " << Nx*Npt << endl;

                    //Memory allocation of the grids in the class TMD
                    F = new double*[Nx];
                    for (int j = 0; j < Nx; ++j)
                        F[j] = new double[Npt];

                    X = new double[Nx];
                    Pt = new double[Npt];
                }
                else
                {
                    istringstream iss(line);
                    iss >> X[i_x] >> Pt[i_pt] >> F[i_x][i_pt];
                    //cout << "\n x = " << X[i_x] << ",\t pt = " << Pt[i_pt] << ", \t F = " << F[i_x][i_pt];

                    i_pt++;
                    if (i_pt == Npt)
                    {
                        i_pt=0;
                        i_x++;
                        //if (i_x == Nx)
                        //    i_x=0;
                    }
                }
            }
        }
    }
    tmd_file.close();

    //Variable range :
    x_min = X[0];
    x_max = X[Nx-1];
    pt_min = Pt[0];
    pt_max = Pt[Npt-1];

    cout << "#\n# Variable range : " << endl;
    cout << "# x_min = " << x_min << ",\tx_max = " << x_max << endl;
    cout << "# pt_min = " << pt_min << ",\tpt_max = " << pt_max << endl;
}

//Create TMD from file
//No header in the file but the number of point in the grid must be given
void MyTMD2::create(string filename, int nx, int npt)
{
    Nx = nx;
    Npt = npt;

    ifstream tmd_file;
    tmd_file.open(filename);

    string line;
    int i_x = 0, i_pt = 0;

    //Memory allocation of the grids in the class TMD
    F = new double*[Nx];
    for (int j = 0; j < Nx; ++j)
        F[j] = new double[Npt];

    X = new double[Nx];
    Pt = new double[Npt];

    if (tmd_file.is_open())
    {
        while (getline(tmd_file, line))
        {
            //if (!(line[0]=='/' && line[1]=='/'))
            //{
                istringstream iss(line);
                iss >> X[i_x] >> Pt[i_pt] >> F[i_x][i_pt];
                //cout << "\n x = " << X[i_x] << ",\t pt = " << Pt[i_pt] << ", \t F = " << F[i_x][i_pt];

                i_pt++;
                if (i_pt == Npt)
                {
                    i_pt=0;
                    i_x++;
                    //if (i_x == Nx)
                    //    i_x=0;
                }
            //}
        }
    }
    tmd_file.close();

    //Variable range :
    x_min = X[0];
    x_max = X[Nx-1];
    pt_min = Pt[0];
    pt_max = Pt[Npt-1];
}

//Create TMD from file
//Here, the file is read once to find Nx and Npt from its structure (as in the PDF routine)
void MyTMD2::create_readN(string filename)
{
    cout << "\n\n#//  Loading TMD from file : " << endl;
    cout << "#///////////////////////////" << endl;
    cout << "# Grid : " << filename << endl;

    Nx = 0;
    Npt = 0;

    findgridN(filename, &Nx, &Npt); //Determines Nx and Npt from file structure

    cout << "#\n# Number of points found : " << endl;
    cout << "# Nx = " << Nx << ", \tNpt = " << Npt << ", \tTotal Nx*Npt = " << Nx*Npt << endl;

    create(filename, Nx, Npt); //Copy grid data in the appropriate class variables

    cout << "#\n# Variable range : " << endl;
    cout << "# x_min = " << x_min << ",\tx_max = " << x_max << endl;
    cout << "# pt_min = " << pt_min << ",\tpt_max = " << pt_max << endl;
}


//Create FF from file
void MyFF::create(string filename)
{
    cout << "\n\n#//  Loading FF from file : " << endl;
    cout << "#//////////////////////////" << endl;
    cout << "# Grid : " << filename << endl;

    ifstream ff_file;
    ff_file.open(filename);

    string line;

    int i_x = 0, i_t = 0;

    bool isNx = false, isNt = false;

    if (ff_file.is_open())
    {
        while (getline(ff_file, line))
        {
            if (line[0]=='/' && line[1]=='/' && line[2]=='N' && line[3]=='x' && line[4]==' ' && line[5]=='=' && line[6]==' ')
            {
                istringstream iss(line);
                string paramname, sign, rest;
                iss >> paramname >> sign >> Nx >> rest;

                isNx = true;

                if (isNt)
                {
                    cout << "#\n# Number of points found : " << endl;
                    cout << "# Nx = " << Nx << ", \tNt = " << Nt << ", \tTotal Nx*Nt = " << Nx*Nt << endl;

                    //Memory allocation of the grids in the class FF
                    D = new double*[Nx];
                    for (int k = 0; k < Nx; ++k)
                        D[k] = new double[Nt];

                    X = new double[Nx];
                    T = new double[Nt];
                }
            }
            else if (line[0]=='/' && line[1]=='/' && line[2]=='N' && line[3]=='t' && line[4]==' ' && line[5]=='=' && line[6]==' ')
            {
                istringstream iss(line);
                string paramname, sign, rest;
                iss >> paramname >> sign >> Nt >> rest;

                isNt = true;

                if (isNx)
                {
                    cout << "#\n# Number of points found : " << endl;
                    cout << "# Nx = " << Nx << ", \tNt = " << Nt << ", \tTotal Nx*Nt = " << Nx*Nt << endl;

                    //Memory allocation of the grids in the class FF
                    D = new double*[Nx];
                    for (int k = 0; k < Nx; ++k)
                        D[k] = new double[Nt];

                    X = new double[Nx];
                    T = new double[Nt];
                }
            }
            else if (line[0]=='/' && line[1]=='/' && line[2]=='e' && line[3]=='p' && line[3]=='s')
            {
                istringstream iss(line);
                string paramname, sign, rest;
                iss >> paramname >> sign >> eps >> rest;
            }
            else if (!(line[0]=='/' && line[1]=='/'))
            {
                if (!isNx || !isNt)
                {
                    cout << "\nError : Nx and Nt not found in header" << endl;
                    exit(-1);
                }
                istringstream iss(line);
                iss >> T[i_t] >> X[i_x] >> D[i_x][i_t];
                //cout << "\n t = " << T[i_t] << ",\t x = " << X[i_x] << ", \t D = " << D[i_x][i_t];

                i_x++;
                if (i_x == Nx)
                {
                    i_x=0;
                    i_t++;
                    //if (i_t == Nt)
                    //    i_t=0;
                }
            }
        }
    }
    ff_file.close();

    //Variable range :
    x_min = X[Nx-1]; // FF grids from Cheb method have x values decreasing
    x_max = X[0];
    t_min = T[0];
    t_max = T[Nt-1];

    cout << "#\n# Variable range : " << endl;
    cout << "# x_min = " << x_min << ",\tx_max = " << x_max << endl;
    cout << "# t_min = " << t_min << ",\tt_max = " << t_max << endl;
}

//Create FF from file
//No header in the file but the number of point in the grid must be given
void MyFF::create(string filename, int nx, int nt)
{
    Nx = nx;
    Nt = nt;

    ifstream ff_file;
    ff_file.open(filename);

    string line;
    int i_x = 0, i_t = 0;

    //Memory allocation of the grids in the class FF
    D = new double*[Nx];
    for (int k = 0; k < Nx; ++k)
        D[k] = new double[Nt];

    X = new double[Nx];
    T = new double[Nt];

    if (ff_file.is_open())
    {
        while (getline(ff_file, line))
        {
            //if (!(line[0]=='/' && line[1]=='/'))
            //{
                istringstream iss(line);
                iss >> T[i_t] >> X[i_x] >> D[i_x][i_t];
                //cout << "\n t = " << T[i_t] << ",\t x = " << X[i_x] << ", \t D = " << D[i_x][i_t];

                i_x++;
                if (i_x == Nx)
                {
                    i_x=0;
                    i_t++;
                    //if (i_t == Nt)
                    //    i_t=0;
                }
            //}
        }
    }
    ff_file.close();

    //Variable range :
    x_min = X[Nx-1]; // FF grids from Cheb method have x values decreasing
    x_max = X[0];
    t_min = T[0];
    t_max = T[Nt-1];
}

//Create FF from file
//Here, the file is read once to find Nx and Npt from its structure (as in the PDF routine)
void MyFF::create_readN(string filename)
{
    cout << "\n\n#//  Loading FF from file : " << endl;
    cout << "#//////////////////////////" << endl;
    cout << "# Grid : " << filename << endl;

    Nx = 0;
    Nt = 0;

    findgridN(filename, &Nt, &Nx); //Determines Nx and Nt from file structure

    cout << "#\n# Number of points found : " << endl;
    cout << "# Nx = " << Nx << ", \tNt = " << Nt << ", \tTotal Nx*Nt = " << Nx*Nt << endl;

    create(filename, Nx, Nt); //Copy grid data in the appropriate class variables

    cout << "#\n# Variable range : " << endl;
    cout << "# x_min = " << x_min << ",\tx_max = " << x_max << endl;
    cout << "# t_min = " << t_min << ",\tt_max = " << t_max << endl;
}



//////////////////////////////////////
// Interpolation related functions  //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Interpolate the value of the PDF distribution on (x, mu2) from its grid (linear interpolation)
double MyPDF::interp(double x, double mu2, parameters param)
{
    int i0, i1, j0, j1;
    double x0, x1, m0, m1;
    double f00, f01, f10, f11;
    //The index / values defined here verify :
    //x0 = X[i0] <= x < x1 = X[i1]
    //m0 = Mu2[j0] <= mu2 < m1 = Mu2[j1]
    //f00 = f[i0][j0] = f(x0, t0)
    //f01 = f[i0][j1] = f(x0, t1)
    //f10 = f[i1][j0] = f(x1, t0)
    //f11 = f[i1][j1] = f(x1, t1)

    //Coefficients of the interpolation
    double X0, X1, M0, M1;

    double fr = -1; // value returned (initiated at unphysical value)

    // out of grid or unphysical values are set to 0
    if ((x > 2*x_max-X[Nx-2]) || (x > 1))         //
        return 0;                                 //
    if ((x < 2*x_min-X[1]) || (x < 0))            //
        return 0;                                 //
    if (mu2 > 2*mu2_max-Mu2[Nmu2-2])              //
        return 0;                                 //
    if ((mu2 < 2*mu2_min-Mu2[1]) || (mu2 < 0))    //
        return 0;                                 //
    ////////////////////////////////////////////////


    // Interpolation in x
    /////////////////////
    if (x > x_max) // values slightly out of the grid are still interpolated
    {
        i0 = Nx-2;
        i1 = Nx-1;
        x0 = X[i0];
        x1 = X[i1];
    }
    else if (x < x_min) // values slightly out of the grid are still interpolated
    {
        i0 = 0;
        i1 = 1;
        x0 = X[i0];
        x1 = X[i1];
    }
    else
    {
        i0 = 0;
        i1 = Nx;
        findinterpindex(&i0, &i1, x, X);
        x0 = X[i0];
        x1 = X[i1];
        //cout << " X[" << i0 << "] = " << x0 << " < x = " << x << " < X[" << i1 << "] = " << x1 << endl;
    }


    // Interpolation in mu2
    ///////////////////////
    if (mu2 > mu2_max) // values slightly out of the grid are still interpolated
    {
        j0 = Nmu2-2;
        j1 = Nmu2-1;
        m0 = Mu2[j0];
        m1 = Mu2[j1];
    }
    else if (mu2 < mu2_min) // values slightly out of the grid are still interpolated
    {
        j0 = 0;
        j1 = 1;
        m0 = Mu2[j0];
        m1 = Mu2[j1];
    }
    else
    {
        j0 = 0;
        j1 = Nmu2;
        findinterpindex(&j0, &j1, mu2, Mu2);
        m0 = Mu2[j0];
        m1 = Mu2[j1];
        //cout << " Mu2[" << j0 << "] = " << m0 << " < mu2 = " << mu2 << " < Mu2[" << j1 << "] = " << m1 << endl;
    }

    f00 = f[i0][j0];
    f01 = f[i0][j1];
    f10 = f[i1][j0];
    f11 = f[i1][j1];

    //Linear interpolation
    //////////////////////
    if (param.interp == 0)
    {
        X0 = (x1-x)/(x1-x0);
        X1 = (x-x0)/(x1-x0);

        M0 = (m1-mu2)/(m1-m0);
        M1 = (mu2-m0)/(m1-m0);

        fr = X0*M0*f00 + X0*M1*f01 + X1*M0*f10 + X1*M1*f11;
        if ((fr < 0)&&(param.positivity)){fr = 0;}
    }

    //Logarithmic interpolation
    ///////////////////////////
    else if (param.interp == 1)
    {
        X0 = (log(x1)-log(x))/(log(x1)-log(x0));
        X1 = (log(x)-log(x0))/(log(x1)-log(x0));

        M0 = (log(m1)-log(mu2))/(log(m1)-log(m0));
        M1 = (log(mu2)-log(m0))/(log(m1)-log(m0));

        fr = X0*M0*f00 + X0*M1*f01 + X1*M0*f10 + X1*M1*f11;
        if ((fr < 0)&&(param.positivity)){fr = 0;}
    }

    //cout << "\nf(x=" << x0 << ",mu2=" << m0 << ") = " << f00 << " ; f(x=" << x << ",mu2=" << mu2 << ") = " << fr << " ; f(x=" << x1 << ",mu2=" << m0 << ") = " << f10 << endl;
    //cout << "f(x=" << x0 << ",mu2=" << m1 << ") = " << f01 << " ; f(x=" << x << ",mu2=" << mu2 << ") = " << fr << " ; f(x=" << x1 << ",mu2=" << m1 << ") = " << f11 << endl;
    return fr;
}

//Interpolate the value of the TMD distribution on (x, pt, mu2) from its grid (linear interpolation)
double MyTMD::interp(double x, double pt, double mu2, parameters param)
{
    int i0, i1, j0, j1, k0, k1;
    double x0, x1, p0, p1, m0, m1;
    double F000, F001, F010, F100, F011, F101, F110, F111;
    //The index / values defined here verify :
    //x0 = X[i0] <= x < x1 = X[i1]
    //p0 = PT[i0] <= pt < p1 = PT[i1]
    //m0 = Mu2[j0] <= mu2 < m1 = Mu2[j1]
    //Fabc = F[ia][jb][kc] = F(xa, pb, mc)

    //Coefficients of the interpolation
    double X0, X1, P0, P1, M0, M1;

    double Fr = -1; // value returned (initiated at unphysical value)

    // out of grid or unphysical values are set to 0
    if ((x > 2*x_max-X[Nx-2]) || (x > 1))         //
        return 0;                                 //
    if ((x < 2*x_min-X[1]) || (x < 0))            //
        return 0;                                 //
    if (pt > 2*pt_max-Pt[Npt-2])                  //
        return 0;                                 //
    if ((pt < 2*pt_max-Pt[1]) || (pt < 0))        //
        return 0;                                 //
    if (mu2 > 2*mu2_max-Mu2[Nmu2-2])              //
        return 0;                                 //
    if ((mu2 < 2*mu2_min-Mu2[1]) || (mu2 < 0))    //
        return 0;                                 //
    ////////////////////////////////////////////////


    // Interpolation in x
    /////////////////////
    if (x > x_max) // values slightly out of the grid are still interpolated
    {
        i0 = Nx-2;
        i1 = Nx-1;
        x0 = X[i0];
        x1 = X[i1];
    }
    else if (x < x_min) // values slightly out of the grid are still interpolated
    {
        i0 = 0;
        i1 = 1;
        x0 = X[i0];
        x1 = X[i1];
    }
    else
    {
        i0 = 0;
        i1 = Nx;
        findinterpindex(&i0, &i1, x, X);
        x0 = X[i0];
        x1 = X[i1];
        //cout << " X[" << i0 << "] = " << x0 << " < x = " << x << " < X[" << i1 << "] = " << x1 << endl;
    }


    // Interpolation in pt
    //////////////////////
    if (pt > pt_max) // values slightly out of the grid are still interpolated
    {
        j0 = Npt-2;
        j1 = Npt-1;
        p0 = Pt[j0];
        p1 = Pt[j1];
    }
    else if (pt < pt_min) // values slightly out of the grid are still interpolated
    {
        j0 = 0;
        j1 = 1;
        p0 = Pt[j0];
        p1 = Pt[j1];
    }
    else
    {
        j0 = 0;
        j1 = Npt;
        findinterpindex(&j0, &j1, pt, Pt);
        p0 = Pt[j0];
        p1 = Pt[j1];
        //cout << " Pt[" << j0 << "] = " << p0 << " < pt = " << pt << " < Pt[" << j1 << "] = " << p1 << endl;
    }


    // Interpolation in mu2
    ///////////////////////
    if (mu2 > mu2_max) // values slightly out of the grid are still interpolated
    {
        k0 = Nmu2-2;
        k1 = Nmu2-1;
        m0 = Mu2[k0];
        m1 = Mu2[k1];
    }
    else if (mu2 < mu2_min) // values slightly out of the grid are still interpolated
    {
        k0 = 0;
        k1 = 1;
        m0 = Mu2[k0];
        m1 = Mu2[k1];
    }
    else
    {
        k0 = 0;
        k1 = Nmu2;
        findinterpindex(&k0, &k1, mu2, Mu2);
        m0 = Mu2[k0];
        m1 = Mu2[k1];
        //cout << " Mu2[" << k0 << "] = " << m0 << " < mu2 = " << mu2 << " < Mu2[" << k1 << "] = " << m1 << endl;
    }


    F000 = F[i0][j0][k0];
    F001 = F[i0][j0][k1];
    F010 = F[i0][j1][k0];
    F100 = F[i1][j0][k0];
    F011 = F[i0][j1][k1];
    F101 = F[i1][j0][k1];
    F110 = F[i1][j1][k0];
    F111 = F[i1][j1][k1];

    //Linear interpolation
    //////////////////////
    if (param.interp == 0)
    {
        X0 = (x1-x)/(x1-x0);
        X1 = (x-x0)/(x1-x0);

        P0 = (p1-pt)/(p1-p0);
        P1 = (pt-p0)/(p1-p0);

        M0 = (m1-mu2)/(m1-m0);
        M1 = (mu2-m0)/(m1-m0);

        Fr = X0*P0*M0*F000 + X1*P0*M0*F100 + X0*P1*M0*F010 + X0*P0*M1*F001 + X1*P1*M0*F110 + X1*P0*M1*F101 + X0*P1*M1*F011 + X1*P1*M1*F111;
        if ((Fr < 0)&&(param.positivity)){Fr = 0;}
    }

    //Logarithmic interpolation
    ///////////////////////////
    else if (param.interp == 1)
    {
        X0 = (log(x1)-log(x))/(log(x1)-log(x0));
        X1 = (log(x)-log(x0))/(log(x1)-log(x0));

        P0 = (log(p1)-log(pt))/(log(p1)-log(p0));
        P1 = (log(pt)-log(p0))/(log(p1)-log(p0));

        M0 = (log(m1)-log(mu2))/(log(m1)-log(m0));
        M1 = (log(mu2)-log(m0))/(log(m1)-log(m0));

        Fr = X0*P0*M0*F000 + X1*P0*M0*F100 + X0*P1*M0*F010 + X0*P0*M1*F001 + X1*P1*M0*F110 + X1*P0*M1*F101 + X0*P1*M1*F011 + X1*P1*M1*F111;
        if ((Fr < 0)&&(param.positivity)){Fr = 0;}
    }

    //cout << "\nf(x=" << x0 << ",pt=" << p0 << ",mu2=" << m0 << ") = " << F000 << " ; f(x=" << x << ",pt=" << pt << ",mu2=" << mu2 << ") = " << Fr << " ; f(x=" << x1 << ",pt=" << p0 << ",mu2=" << m0 << ") = " << F100 << endl;
    //cout << "f(x=" << x0 << ",pt=" << p0 << ",mu2=" << m1 << ") = " << F001 << " ; f(x=" << x << ",pt=" << pt << ",mu2=" << mu2 << ") = " << Fr << " ; f(x=" << x1 << ",pt=" << p0 << ",mu2=" << m1 << ") = " << F101 << endl;
    //cout << "f(x=" << x0 << ",pt=" << p1 << ",mu2=" << m1 << ") = " << F011 << " ; f(x=" << x << ",pt=" << pt << ",mu2=" << mu2 << ") = " << Fr << " ; f(x=" << x1 << ",pt=" << p1 << ",mu2=" << m1 << ") = " << F111 << endl;
    return Fr;
}

//Interpolate the value of the TMD distribution on (x, pt, mu2) from its grid (linear interpolation)
double MyTMD2::interp(double x, double pt, parameters param)
{
    int i0, i1, j0, j1;
    double x0, x1, p0, p1;
    double F00, F01, F10, F11;
    //The index / values defined here verify :
    //x0 = X[i0] <= x < x1 = X[i1]
    //p0 = PT[i0] <= pt < p1 = PT[i1]
    //Fabc = F[ia][jb] = F(xa, pb)

    //Coefficients of the interpolation
    double X0, X1, P0, P1;

    double Fr = -1; // value returned (initiated at unphysical value)

    //cout << "x = " << x << ",\t2*x_max-X[Nx-2] = " << 2*x_max-X[Nx-2] << ",\t2*x_min-X[1] = " << 2*x_min-X[1] << endl;
    //cout << "pt = " << pt << ",\t2*pt_max-Pt[Npt-2] = " << 2*pt_max-Pt[Npt-2] << ",\t2*pt_min-Pt[1] = " << 2*pt_min-Pt[1] << endl;

    // out of grid or unphysical values are set to 0
    ///////////////////////////////////////////////////////////
    if ((x > 2*x_max-X[Nx-2]) || (x < 2*x_min-X[1]))         //
        return 0;                                            //
    if ((x < 0) || (x > 1))                                  //
        return 0;                                            //
    if ((pt < 2*pt_min-Pt[1]) || (pt > 2*pt_max-Pt[Npt-2]))  //
        return 0;                                            //
    if (pt < 0)                                              //
        return 0;                                            //
    ///////////////////////////////////////////////////////////


    // Interpolation in x
    /////////////////////
    if (x > x_max) // values slightly out of the grid are still interpolated
    {
        i0 = Nx-2;
        i1 = Nx-1;
        x0 = X[i0];
        x1 = X[i1];
    }
    else if (x < x_min) // values slightly out of the grid are still interpolated
    {
        i0 = 0;
        i1 = 1;
        x0 = X[i0];
        x1 = X[i1];
    }
    else
    {
        i0 = 0;
        i1 = Nx;
        findinterpindex(&i0, &i1, x, X);
        x0 = X[i0];
        x1 = X[i1];
        //cout << " X[" << i0 << "] = " << x0 << " < x = " << x << " < X[" << i1 << "] = " << x1 << endl;
    }


    // Interpolation in pt
    //////////////////////
    if (pt > pt_max) // values slightly out of the grid are still interpolated
    {
        j0 = Npt-2;
        j1 = Npt-1;
        p0 = Pt[j0];
        p1 = Pt[j1];
    }
    else if (pt < pt_min) // values slightly out of the grid are still interpolated
    {
        j0 = 0;
        j1 = 1;
        p0 = Pt[j0];
        p1 = Pt[j1];
    }
    else
    {
        j0 = 0;
        j1 = Npt;
        findinterpindex(&j0, &j1, pt, Pt);
        p0 = Pt[j0];
        p1 = Pt[j1];
        //cout << " Pt[" << j0 << "] = " << p0 << " < pt = " << pt << " < Pt[" << j1 << "] = " << p1 << endl;
    }


    F00 = F[i0][j0];
    F01 = F[i0][j1];
    F10 = F[i1][j0];
    F11 = F[i1][j1];

    //Linear interpolation
    //////////////////////
    if (param.interp == 0)
    {
        X0 = (x1-x)/(x1-x0);
        X1 = (x-x0)/(x1-x0);

        P0 = (p1-pt)/(p1-p0);
        P1 = (pt-p0)/(p1-p0);

        Fr = X0*P0*F00 + X1*P0*F10 + X0*P1*F01 + X1*P1*F11;
        if ((Fr < 0)&&(param.positivity)){Fr = 0;}
    }

    //Logarithmic interpolation
    ///////////////////////////
    else if (param.interp == 1)
    {
        X0 = (log(x1)-log(x))/(log(x1)-log(x0));
        X1 = (log(x)-log(x0))/(log(x1)-log(x0));

        P0 = (log(p1)-log(pt))/(log(p1)-log(p0));
        P1 = (log(pt)-log(p0))/(log(p1)-log(p0));

        Fr = X0*P0*F00 + X1*P0*F10 + X0*P1*F01 + X1*P1*F11;
        if ((Fr < 0)&&(param.positivity)){Fr = 0;}
    }

    //cout << "\nF(x=" << x0 << ",pt=" << p0 << ") = " << F00 << " ; F(x=" << x << ",pt=" << pt << ") = " << Fr << " ; F(x=" << x1 << ",pt=" << p0 << ") = " << F10 << endl;
    //cout << "F(x=" << x0 << ",pt=" << p1 << ") = " << F01 << " ; F(x=" << x << ",pt=" << pt << ") = " << Fr << " ; F(x=" << x1 << ",pt=" << p1 << ") = " << F11 << endl;
    return Fr;
}

//Interpolate the value of the FF distribution on (x, t) from its grid (linear interpolation)
double MyFF::interp(double x, double t, parameters param)
{
    int i0, i1, j0, j1;
    double x0, x1, t0, t1;
    double D00, D01, D10, D11;
    //The index / values defined here verify :
    //x0 = X[i0] <= x < x1 = X[i1]
    //t0 = T[j0] <= t < t1 = T[j1]
    //D00 = D[i0][j0] = D(x0, t0)
    //D01 = D[i0][j1] = D(x0, t1)
    //D10 = D[i1][j0] = D(x1, t0)
    //D11 = D[i1][j1] = D(x1, t1)

    //Coefficients of the interpolation
    double X0, X1, T0, T1;

    double Dr = -1; // value returned (initiated at unphysical value)

    // out of grid or unphysical values are set to 0
    ////////////////////////////////////////////////////
    if ((x < 2*x_min-X[Nx-2]) || (x > 2*x_max-X[1]))  //
        return 0;                                     //
    if ((x < 0) || (x > 1))                           //
        return 0;                                     //
    if ((t < 2*t_min-T[1]) || (t > 2*t_max-T[Nt-2]))  //
        return 0;                                     //
    if (t < 0)                                        //
        return 0;                                     //
    ////////////////////////////////////////////////////


    // Interpolation in x
    /////////////////////
    if (x > x_max) // values slightly out of the grid are still interpolated
    {
        i0 = 0;
        i1 = 1;
        x0 = X[i0];
        x1 = X[i1];
    }
    else if (x < x_min) // values slightly out of the grid are still interpolated
    {
        i0 = Nx-2;
        i1 = Nx-1;
        x0 = X[i0];
        x1 = X[i1];
    }
    else
    {
        i0 = 0;
        i1 = Nx;
        findinterpindex_inv(&i0, &i1, x, X);
        x0 = X[i0];
        x1 = X[i1];
        //cout << " X[" << i0 << "] = " << x0 << " < x = " << x << " < X[" << i1 << "] = " << x1 << endl;
    }


    // Interpolation in t
    /////////////////////
    if (t > t_max) // values slightly out of the grid are still interpolated
    {
        j0 = Nt-2;
        j1 = Nt-1;
        t0 = T[j0];
        t1 = T[j1];
    }
    else if (t < t_min) // values slightly out of the grid are still interpolated
    {
        j0 = 0;
        j1 = 1;
        t0 = T[j0];
        t1 = T[j1];
    }
    else
    {
        j0 = 0;
        j1 = Nt;
        findinterpindex(&j0, &j1, t, T);
        t0 = T[j0];
        t1 = T[j1];
        //cout << " T[" << j0 << "] = " << t0 << " < t = " << t << " < T[" << j1 << "] = " << t1 << endl;
    }


    D00 = D[i0][j0];
    D01 = D[i0][j1];
    D10 = D[i1][j0];
    D11 = D[i1][j1];

    //Linear interpolation
    //////////////////////
    if (param.interp == 0)
    {
        X0 = (x1-x)/(x1-x0);
        X1 = (x-x0)/(x1-x0);

        T0 = (t1-t)/(t1-t0);
        T1 = (t-t0)/(t1-t0);

        Dr = X0*T0*D00 + X0*T1*D01 + X1*T0*D10 + X1*T1*D11;
        if ((Dr < 0)&&(param.positivity)){Dr = 0;}
    }

    //Logarithmic interpolation
    ///////////////////////////
    else if (param.interp == 1)
    {
        X0 = (log(x1)-log(x))/(log(x1)-log(x0));
        X1 = (log(x)-log(x0))/(log(x1)-log(x0));

        T0 = (log(t1)-log(t))/(log(t1)-log(t0));
        T1 = (log(t)-log(t0))/(log(t1)-log(t0));

        Dr = X0*T0*D00 + X0*T1*D01 + X1*T0*D10 + X1*T1*D11;
        if ((Dr < 0)&&(param.positivity)){Dr = 0;}
    }

    //cout << "\nD(x=" << x0 << ",t=" << t0 << ") = " << D00 << " ; D(x=" << x << ",t=" << t << ") = " << Dr << " ; D(x=" << x1 << ",t=" << t0 << ") = " << D10 << endl;
    //cout << "D(x=" << x0 << ",t=" << t1 << ") = " << D01 << " ; D(x=" << x << ",t=" << t << ") = " << Dr << " ; D(x=" << x1 << ",t=" << t1 << ") = " << D11 << endl;
    return Dr;
}



/////////////////////////////////////////
// Additional functions / subroutines  //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Recursive function to find i0 and i1 such as X[i0] < x < X[i1], i1 = i0+1 by dychotomy, giving i0 = 0, i1 = Nx
void findinterpindex(int *i0, int *i1, double x, double *&X)
{
    //cout << "\nInterpolation : X[i0=" << *i0 << "]=" << X[*i0] << " < x=" << x << " < X[i1=" << *i1 << "]=" << X[*i1];

    if (*i1 == *i0 + 1)
        return;
    else if (x <= X[*i0 + (*i1-*i0)/2])
    {
        *i1 = *i0 + (*i1-*i0)/2;
        findinterpindex(i0, i1, x, *&X);
    }
    else
    {
        *i0 = *i0 + (*i1-*i0)/2;
        findinterpindex(i0, i1, x, *&X);
    }
}

// Recursive function to find i0 and i1 such as X[i1] < x < X[i0], i0 = i1+1 by dychotomy, giving i0 = Nx, i1 = 0
void findinterpindex_inv(int *i0, int *i1, double x, double *&X)
{
    //cout << "\nInterpolation : X[i0=" << *i0 << "]=" << X[*i0] << " > x=" << x << " > X[i1=" << *i1 << "]=" << X[*i1];

    if (*i1 == *i0 + 1)
        return;
    else if (x < X[*i0 + (*i1-*i0)/2])
    {
        *i0 = *i0 + (*i1-*i0)/2;
        findinterpindex_inv(i0, i1, x, *&X);
    }
    else
    {
        *i1 = *i0 + (*i1-*i0)/2;
        findinterpindex_inv(i0, i1, x, *&X);
    }
}

// find the number of point x and y in a grid of the form : x y f(x,y)
// It is supposed that the grid is organized with y increasing for fixed values of x then x increasing (while y does Nx time the same loop)
void findgridN(string filename, int *Nx, int *Ny)
{
    ifstream grid_file;
    grid_file.open(filename);

    string line;
    int init = 0, i_x = 1, i_y = 1;

    double x[2], y[2], D;

    if (grid_file.is_open())
    {
        while (getline(grid_file, line))
        {
            if (!(line[0]=='/' && line[1]=='/'))
            {
                istringstream iss(line);
                iss >> x[1] >> y[1] >> D;
                //cout << x[1] << "\t" << y[1] << "\t" << D << endl;

                if (init == 1)
                {
                    if (!(y[0] == y[1]) && !(x[0] == x[1]))
                    {
                        i_x += 1;
                        i_y = 1;
                    }
                    else if (!(y[0] == y[1]))
                        i_y += 1;
                }

                x[0] = x[1];
                y[0] = y[1];
                init = 1;
            }
        }
    }
    *Nx = i_x;
    *Ny = i_y;

    //cout << "Number of points found : " << endl;
    //cout << " Nx = " << *Nx << ", \tNy = " << *Ny << ", \tTotal Nx.Ny = " << *Ny * *Nx << endl;

    return;
}

// find the number of point x, and z y in a grid of the form : x y z f(x,y,z)
// It is supposed that the grid is organized with y increasing for fixed values of x then x increasing (while y does Nx time the same loop)
void findgridN(string filename, int *Nx, int *Ny, int *Nz)
{
    ifstream ff_file;
    ff_file.open(filename);

    string line;
    int init = 0, i_x = 1, i_y = 1, i_z = 1;

    double x[2], y[2], z[2], D;

    if (ff_file.is_open())
    {
        while (getline(ff_file, line))
        {
            if (!(line[0]=='/' && line[1]=='/'))
            {
                istringstream iss(line);
                iss >> x[1] >> y[1] >> z[1] >> D;
                //cout << x[1] << "\t" << y[1] << "\t" << z[1] << "\t" << D << endl;

                if (init == 1)
                {
                    if (!(z[0] == z[1]) && !(y[0] == y[1]) && !(x[0] == x[1]))
                    {
                        i_x += 1;
                        i_y = 1;
                        i_z = 1;
                    }
                    else if (!(z[0] == z[1]) && !(y[0] == y[1]))
                    {
                        i_y += 1;
                        i_z = 1;
                    }
                    else if (!(z[0] == z[1]))
                        i_z += 1;
                }

                x[0] = x[1];
                y[0] = y[1];
                z[0] = z[1];
                init = 1;
            }
        }
    }
    *Nx = i_x;
    *Ny = i_y;
    *Nz = i_z;

    //cout << "Number of points found : " << endl;
    //cout << "\tNx = " << *Nx << ", \tNy = " << *Ny << ", \tNz = " << *Nz << ", \tTotal Nx*Ny*Nz = " << *Nx * *Ny * *Nz << endl;

    return;
}



/////////////////////////////////
// Destructors of the classes  //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
//PDF destructor
MyPDF::~MyPDF()
{
    for (int i = 0; i < Nx; ++i)
        delete[] f[i];
    delete[] f;
    delete[] X;
    delete[] Mu2;
}

//TMD destructor
MyTMD::~MyTMD()
{
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Npt; ++j)
            delete[] F[i][j];
        delete[] F[i];
    }
    delete[] F;
    delete[] X;
    delete[] Pt;
    delete[] Mu2;
}

//TMD destructor
MyTMD2::~MyTMD2()
{
    for (int i = 0; i < Nx; ++i)
        delete[] F[i];
    delete[] F;
    delete[] X;
    delete[] Pt;
}

//FF destructor
MyFF::~MyFF()
{
    for (int i = 0; i < Nx; ++i)
        delete[] D[i];
    delete[] D;
    delete[] X;
    delete[] T;
}//*/

#endif // DISTRIBUTION_CLASSES_H
