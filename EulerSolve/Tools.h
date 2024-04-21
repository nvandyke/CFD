#ifndef TOOLS_H
#define TOOLS_H

#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "matrix.h"
#include "block.h"

using namespace std;
const double y = 1.4;
const double pi = 3.14159265358979323846;

enum BoundaryCondition {
    Inviscid_Wall,
    Subsonic_Outlet,
    Subsonic_Inlet,
    Supersonic_Outlet,
    Supersonic_Inlet,
    Freestream
};

enum Flux {
    Roe,
    HLLE,
    Rusanov
};

void printResults(Matrix& u, Matrix& e);

Matrix roots(double a, double b, double c);

double max(double a, double b);
double min(double a, double b);

double waveSpeed(Matrix& u, Matrix& n);

Matrix flux(Matrix& uL, Matrix& uR, Matrix& n, Flux R);

Matrix F(Matrix& u);

Matrix Outflow(Matrix& u, Matrix& n, double Pb);

Matrix Inflow(Matrix& u, Matrix& n, double Tt, double Pt, double a, double R);

Matrix wallFlux(Matrix& u, Matrix& n);


//Finite Volume Storage
class FVConditions {
public:
    double M, a, R, Pinf, Tt, Pt, c, p;
    Matrix v;
    Matrix uinf;
    void set();
    FVConditions();
    FVConditions(double Mach, double aoa);
};


class FVmesh {
public:
    Block I2E;
    Block B2E;
    Matrix In;
    Matrix Bn;
    Matrix Il;
    Matrix Bl;
    Matrix A;
    Matrix Br;
    Matrix Ir;
    Matrix C;

    Matrix V;
    Block E;
    Block* B;
    string* Bname;

    FVmesh(string filename);
    Block connectivity(Block& E);
    vector<int> edge2vertex(int t, int e, Block& E);
    Block boundaryConnectivity(Block& E, Block& boundary, int bgroup);
    Matrix normal(int t, int e);
    double length(int t, int e);
    Matrix midpoint(int t, int e);
    Matrix pointsFromTE(int t, int e);
    double Area(int t);
    double verify();


};


class FVstate {
public:
    Matrix u;
    int Order;
    Matrix gradux;
    Matrix graduy;
    FVstate();
    FVstate(int O, FVmesh m, FVConditions c);
};


//Finite Volume Methods
Matrix FV_solve(FVstate& u, FVmesh& m, FVConditions c);

Matrix residual(FVstate& u, FVmesh& m, FVConditions& c, Matrix& dt);

void gradient(FVstate& u, FVmesh& m);


//Finite Element Storage
class FEConditions 
{
public:
	double M, a, R, Pinf, Tt, Pt, c, p;
	Matrix v;
	Matrix uinf;




};


class FEmesh
{
public:


};


class FEstate
{
public:
	int p;//numerical order of accuracy
	int q;//geometry order of accuraycy
    Matrix u;

	FEstate() {
		p = 1;
		q = 1;
	};

	FEstate(int pp, int qq, const FEmesh& m, const FEConditions& c) {

	};


};


//Finite Element Methods


#endif