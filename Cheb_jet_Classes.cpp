#include "Classes.h"
#include <stdlib.h>

//Create PDF from file
void PDF::create(string filename)
{
    ifstream pdf_file;
    pdf_file.open(filename);

    string line;
    int i = 0;
    int i_x, i_mu2, i_f;
    double dt;

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

                    //Memory allocation of the grids in the class FF
                    f = new double*[Nx];
                    for (int k = 0; k < Nx; ++k)
                        f[i] = new double[Nmu2];

                    x = new double[Nx];
                    mu2 = new double[Nmu2];
                }
                else
                {
                    istringstream ss(line);
                    ss >> mu2[i_mu2] >> x[i_x] >> f[i_x][i_mu2];
                    cout << "\n mu2 = " << mu2[i_mu2] << ",\t x = " << x[i_x] << ", \t f = " << f[i_x][i_mu2];
                }
            }
        }
    }
}

//Delete PDF
/*
PDF::~PDF()
{

}*/
/*
//Create FF from file
FF::FF(string filename)
{
    ifstream pdf_file;
    pdf_file.open(filename);

    string line;
    int i = 0;
    double dt;

    if (parameters_file.is_open())
    {
        while (getline(parameters_file, line))
        {
            if (!(lineInput[0]=='/' && lineInput[1]=='/'))
            {
                if (i==0)
                {
                    Nx = strtod(line.c_str(), NULL);
                    i++;
                }
                else if (i==1)
                {
                    Nt = strtod(line.c_str(), NULL);
                    i++;
                }
                else if (i==2)
                {
                    dt = strtod(line.c_str(), NULL);
                    cout << "\n-n_t = " << *n_t;
                    i++;

                    //Memory allocation of the grids in the class FF
                    D = new double*[Nx];
                    for (int k = 0; k < Nx; ++k)
                        D[i] = new double[Nt];

                    x = new double*[Nx];
                    t = new double*[Nt];
                }
                else
                {

                }
            }
        }
    }
}*/
