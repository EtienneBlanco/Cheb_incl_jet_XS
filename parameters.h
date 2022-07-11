/*=====================================================\
|                                                      |
|                      Parameters                      |
|	                                                   |
\=====================================================*/


// On this file, a parameters class is defined with some default values and routines to read them from a file.

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <sstream>

//const long double M_PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
double paramstr2double(string str); //Function that find a double in a parameter line and return it
int paramstr2int(string str); //Function that find an integer in a parameter and return it
bool paramstr2bool(string str); //Function that find a boolean in a parameter and return it
string paramstr2str(string str); //Function that find a string in a parameter and return it

//Constants parameters used
class parameters
{
    public:
        int evo; //type of evolution : pure gluon (0) or with quarks (1)
        int intmethod; //integration method (0 for Trapezoidal, 1 for Chebyshev, 2-3 : Left-Right Riemann sum)
        int interp; // Method to use for the interpolation (0: linear, 1: log)
        bool positivity; // Consider or not negative values in distribution grids (if equal to 1, all negative values are set to 0)

        int Ny; //Number of rapidity grid points
        int Np; //Number of transverse momentum grid points
        int Nx; //Number of Chebychev nodes for the integration in x

        double ymin; //Minimum rapidity value
        double ymax; //Maximum rapidity value
        double ptmin; //Minimum transverse momentum
        double ptmax; //Maximum transverse momentum
        double qtmin; //Minimum hard transverse momentum
        double qtmax; //Maximum hard transverse momentum

        double alphaS; //Coupling constant
        double alphabar; // from grid?
        double S; //Center of mass energy of the collision (actually sqrt(s))
        double Temp; //Temperature of the medium
        double qhat; //Quenching parameter (from grid?)
        double t; //lenght of the medium [fm]
        double eps; //low-x cut-off to use logarithmic scale in x
        int Nc; //Number of colors
        double Cf;

        // Internal parameters  (not set to be changed with the 'Parameter' input file)
        bool paramfromfiles;
        int Nf = 3; //Number of quark flavor (for LHAPDF)
        double fmtoGeV = 0.19732; // fermi to GeV-1
        double pbtoGeV = 2.56819*1e-9; // pico barn to GeV-2 conversion

        string PDFset; //name of the PDF set (LHAPDF)
        string TMDsetga; //name of the TMD grid for gluon (adjoint)
        string TMDsetgf; //name of the TMD grid for gluon (fundamental)
        string FFsetg2g; //name of the FF grid for gluons in gluon jet
        string FFsetg2S; //name of the FF grid for quark singlets in gluon jet
        string FFsetS2S; //name of the FF grid for quark singlets in quark singlet jet
        string FFsetS2g; //name of the FF grid for gluons in quark singlet jet

        string gridname; //name to be used for the result grid

        void from_file(string filename); //initialize parameters from a parameter file
        void default_(); //initialize parameters with default values
        void show(string parameterfile); //Display parameters in the class
};


//function reading the parameters file
void parameters::from_file(string filename)
{
    string line;

    ifstream parameters_file;
    parameters_file.open(filename);

    default_(); // Initialize the parameters with default values in case some are not defined in file

    if (parameters_file.is_open())
    {
        while (getline(parameters_file, line))
        {
            if (!(line[0]=='/' && line[1]=='/'))
            {

                if (line[0]=='e' && line[1]=='v' && line[2]=='o')
                {
                    evo = paramstr2int(line);
                    if (!(evo == 0) && !(evo == 1))
                        evo = 1;
                }

                if (line[0]=='i' && line[1]=='n' && line[2]=='t' && line[3]=='e' && line[4]=='g' && line[5]=='r' && line[6]=='a' && line[7]=='t' && line[8]=='i' && line[9]=='o' && line[10]=='n')
                {
                    intmethod = paramstr2int(line);
                    if (!(intmethod == 0) && !(intmethod == 1) && !(intmethod == 2) && !(intmethod == 3))
                        intmethod = 0;
                }
                if (line[0]=='i' && line[1]=='n' && line[2]=='t' && line[3]=='e' && line[4]=='r' && line[5]=='p')
                {
                    interp = paramstr2int(line);
                    if (!(interp == 0) && !(interp == 1))
                        interp = 0;
                }
                if (line[0]=='p' && line[1]=='o' && line[2]=='s' && line[3]=='i' && line[4]=='t' && line[5]=='i' && line[6]=='v' && line[7]=='i' && line[8]=='t' && line[9]=='y')
                    positivity = paramstr2bool(line);

                if (line[0]=='N' && line[1]=='y')
                    Ny = paramstr2int(line);
                if (line[0]=='N' && line[1]=='p')
                    Np = paramstr2int(line);
                if (line[0]=='N' && line[1]=='x')
                    Nx = paramstr2int(line);
                if (line[0]=='e' && line[1]=='p' && line[2]=='s')
                    eps = paramstr2double(line);

                if (line[0]=='y' && line[1]=='m' && line[2]=='i' && line[3]=='n')
                    ymin = paramstr2double(line);
                if (line[0]=='y' && line[1]=='m' && line[2]=='a' && line[3]=='x')
                    ymax = paramstr2double(line);
                if (line[0]=='p' && line[1]=='t' && line[2]=='m' && line[3]=='i' && line[4]=='n')
                    ptmin = paramstr2double(line);
                if (line[0]=='p' && line[1]=='t' && line[2]=='m' && line[3]=='a' && line[4]=='x')
                    ptmax = paramstr2double(line);
                if (line[0]=='q' && line[1]=='t' && line[2]=='m' && line[3]=='i' && line[4]=='n')
                    qtmin = paramstr2double(line);
                if (line[0]=='q' && line[1]=='t' && line[2]=='m' && line[3]=='a' && line[4]=='x')
                    qtmax = paramstr2double(line);

                if (line[0]=='a' && line[1]=='l' && line[2]=='p' && line[3]=='h' && line[4]=='a' && line[5]=='S')
                {
                    alphaS = paramstr2double(line);
                    alphabar = alphaS/M_PI;
                }
                if (line[0]=='a' && line[1]=='l' && line[2]=='p' && line[3]=='h' && line[4]=='a' && line[5]=='b' && line[6]=='a' && line[7]=='r')
                    alphabar = paramstr2double(line);
                if (line[0]=='S')
                    S = paramstr2double(line);
                if (line[0]=='T' && line[1]=='e' && line[2]=='m' && line[3]=='p')
                    Temp = paramstr2double(line);
                if (line[0]=='q' && line[1]=='h' && line[2]=='a' && line[3]=='t')
                    qhat = paramstr2double(line);
                if (line[0]=='t')
                    t = paramstr2double(line);
                if (line[0]=='N' && line[1]=='c')
                {
                    Nc = paramstr2int(line);
                    Cf = (Nc*Nc-1.)/2./Nc;
                }
                if (line[0]=='C' && line[1]=='f')
                    Cf = paramstr2double(line);

                if (line[0]=='P' && line[1]=='D' && line[2]=='F' && line[3]=='s' && line[4]=='e' && line[5]=='t')
                    PDFset = paramstr2str(line);
                if (line[0]=='T' && line[1]=='M' && line[2]=='D' && line[3]=='s' && line[4]=='e' && line[5]=='t' && line[6]=='g' && line[7]=='a')
                    TMDsetga = paramstr2str(line);
                if (line[0]=='T' && line[1]=='M' && line[2]=='D' && line[3]=='s' && line[4]=='e' && line[5]=='t' && line[6]=='g' && line[7]=='f')
                    TMDsetgf = paramstr2str(line);
                if (line[0]=='F' && line[1]=='F' && line[2]=='s' && line[3]=='e' && line[4]=='t' && line[5]=='g' && line[6]=='2' && line[7]=='g')
                    FFsetg2g = paramstr2str(line);
                if (line[0]=='F' && line[1]=='F' && line[2]=='s' && line[3]=='e' && line[4]=='t'&& line[5]=='S' && line[6]=='2' && line[7]=='g')
                    FFsetS2g = paramstr2str(line);
                if (line[0]=='F' && line[1]=='F' && line[2]=='s' && line[3]=='e' && line[4]=='t'&& line[5]=='S' && line[6]=='2' && line[7]=='S')
                    FFsetS2S = paramstr2str(line);
                if (line[0]=='F' && line[1]=='F' && line[2]=='s' && line[3]=='e' && line[4]=='t'&& line[5]=='g' && line[6]=='2' && line[7]=='S')
                    FFsetg2S = paramstr2str(line);
                if (line[0]=='f' && line[1]=='i' && line[2]=='l' && line[3]=='e' && line[4]=='n' && line[5]=='a' && line[6]=='m' && line[7]=='e')
                    gridname = paramstr2str(line);
            }
        }
        paramfromfiles = true;
    }
    else
        paramfromfiles = false;
    show(filename); //display the parameters that will be used

    parameters_file.close();
}

//Default setting for the parameters
void parameters::default_()
{
    evo = 1;
    intmethod = 0;
    interp = 0;
    positivity = 0;

    Ny = 50;
    Np = 50;
    Nx = 50;

    ymin = 0;
    ymax = 4;
    ptmin = 0;
    ptmax = 200;
    qtmin = 1;
    qtmax = 200;

    alphaS = 0.1*M_PI;
    alphabar = alphaS/M_PI;
    S = 5000.; //[GeV]
    Temp = 250; // [MeV]
    qhat = 1; // [GeV^2/fm]
    t = 0.1; // [fm]
    eps = 1e-5;
    Nc = 3;
    Cf = (Nc*Nc-1.)/2./Nc;

    PDFset = "CT10nlo";
    TMDsetga = "../Distributions/TMD/glueadjointKS.dat";
    TMDsetgf = "../Distributions/TMD/glueadfundamentalKS.dat";
    FFsetg2g = "../Distributions/FF/cheb_g2g_x_tau";
    FFsetg2S = "../Distributions/FF/cheb_g2S_x_tau";
    FFsetS2S = "../Distributions/FF/cheb_S2S_x_tau";
    FFsetS2g = "../Distributions/FF/cheb_S2g_x_tau";

    gridname = "xs";
}

//Display the parameters
void parameters::show(string filename)
{
    cout << "\n#//////////////////#" << endl;
    cout << "#//  Parameters  //#" << endl;
    cout << "#/////////////////////////////////////////////#" << endl;
    if (paramfromfiles)
    {
        cout << "# Parameters are read from the following file : " << filename << endl;
        cout << "# and are set to the following values : " << endl;
    }
    else
    {
        cout << "# /!\\ 'Parameters' file not found" << endl;
        cout << "# Parameters are set to the following default values : " << endl;
    }
    cout << "#" << endl;

    cout << "# -Evolution equations : ";
    if (evo == 0)
        cout << "Pure gluons" << endl;
    else if (evo == 1)
        cout << "Quarks and gluons" << endl;

    cout << "# -Integration method : ";
    if (intmethod == 0)
        cout << "Trapezoid" << endl;
    else if (intmethod == 1)
        cout << "Gauss-Chebyshev" << endl;
    else if (intmethod == 2)
        cout << "Left Riemann sum" << endl;
    else if (intmethod == 3)
        cout << "Right Riemann sum" << endl;

    cout << "# -Interpolation method : ";
    if (interp == 0)
        cout << "Linear" << endl;
    else if (interp == 1)
        cout << "Logarithmic" << endl;

    cout << "# -Positivity : ";
    if (positivity)
        cout << "Yes" << endl;
    else
        cout << "No" << endl;

    cout << "#" << endl;

    cout << "# -Ny = " << Ny << endl;
    cout << "# -Npt = " << Np << endl;
    cout << "# -Nx = " << Nx << endl;

    cout << "#" << endl;

    cout << "# -ymin = " << ymin << endl;
    cout << "# -ymax = " << ymax << endl;
    cout << "# -ptmin = " << ptmin << " [GeV]" << endl;
    cout << "# -ptmax = " << ptmax << " [GeV]" << endl;
    cout << "# -qtmin = " << qtmin << " [GeV]" << endl;
    cout << "# -qtmax = " << qtmax << " [GeV]" << endl;

    cout << "#" << endl;

    cout << "# -alphaS = " << alphaS << endl;
    cout << "# -alphabar = " << alphabar << endl;
    cout << "# -S = " << S << endl;
    cout << "# -T = " << Temp << " [MeV]" << endl;
    cout << "# -qhat = " << qhat << " [GeV^2/fm]" << endl;
    cout << "# -t = " << t << " [fm]" << endl;
    cout << "# -eps = " << eps << endl;
    cout << "# -Nc = " << Nc << endl;
    cout << "# -Cf = " << Cf << endl;

    cout << "#" << endl;

    cout << "# -PDF set = " << PDFset << endl;
    cout << "# -TMD set for gluons (adjoint) = " << TMDsetga << endl;
    cout << "# -TMD set for gluons (fundamental) = " << TMDsetgf << endl;
    cout << "# -FF set for gluons in gluon jet = " << FFsetg2g << endl;
    cout << "# -FF set for gluons in S jet = " << FFsetS2g << endl;
    cout << "# -FF set for S in S jet = " << FFsetS2S << endl;
    cout << "# -FF set for S in gluon jet = " << FFsetg2S << endl;
    cout << "# (S stands for quark singlet)" << endl;

    cout << "#" << endl;

    cout << "# -gridname = " << gridname << endl;
    cout << "#/////////////////////////////////////////////#" << endl;
}

double paramstr2double(string str)
{
    istringstream iss(str);
    string paramname, sign, rest;
    double param;

    iss >> paramname >> sign >> param >> rest;
    return param;
}

int paramstr2int(string str)
{
    istringstream iss(str);
    string paramname, sign, rest;
    int param;

    iss >> paramname >> sign >> param >> rest;
    return param;
}

bool paramstr2bool(string str)
{
    istringstream iss(str);
    string paramname, sign, rest;
    bool param;

    iss >> paramname >> sign >> param >> rest;
    return param;
}

string paramstr2str(string str)
{
    istringstream iss(str);
    string paramname, sign, rest;
    string param;

    iss >> paramname >> sign >> param >> rest;
    return param;
}

#endif // CHEB_XS_H
