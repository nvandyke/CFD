#ifndef TOOLS_H
#define TOOLS_H

#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "matrix.h"
#include "block.h"

const double pi = 3.14159265358979323846;
const double y = 1.3;

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

double waveSpeed(Matrix& u, Matrix& n);

Matrix flux(Matrix& uL, Matrix& uR, Matrix& n, int R);

Matrix F(Matrix& u);

Matrix Outflow(Matrix& u, Matrix& n, double Pb);

Matrix Inflow(Matrix& u, Matrix& n, double Tt, double Pt, double a, double R);

Matrix wallFlux(Matrix& u, Matrix& n);


class Mesh {
private:
    int* n_1d;
    double* x_1d;
    double* w_1d;

    int* n_2d;
    double* x_2d;
    double* w_2d;

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

    Mesh();
    Mesh(std::string filename);
    Block connectivity();
    std::vector<int> edge2vertex(int t, int e);
    Block boundaryConnectivity(Block& boundary, int bgroup);
    Matrix normal(int t, int e);
    double length(int t, int e);
    Matrix midpoint(int t, int e);
    Matrix pointsFromTE(int t, int e);
    double Area(int t);
    void Jacobians1d(int elem);
    double* RefEdge2Elem(int edge);
    double verify();
    void quad1d(int order);
    void quad2d(int order);


};


class State {
private:
    double M, a, Rinf, Pinf, Tt, Pt, c, p;
    Matrix v;
    Matrix uinf;
    Mesh mesh;
    void set();
    Matrix gradux;
    Matrix graduy;
    const double CFL = 0.5;
    double tol = 1e-7;
    double max = 1 / tol;
    int MaxIter = 100000;
    int flux_val = Roe;
public:
    Matrix u;
    int Order, Order_geom;
    State();
    State(int O, Mesh m, double Mach, double aoa);

    //Finite Volume Methods
    Matrix FV_solve();

    Matrix residual(Matrix& dt);
    void interioredges(double* R, double* sl);
    void boundaryedges(double* R, double* sl);

    void gradient();
    void grad_internal(double* gradux, double* graduy);
    void grad_boundary(double* gradux, double* graduy);

    //Finite Element Methods
    Matrix FE_solve();
};

#endif