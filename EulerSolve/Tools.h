#ifndef TOOLS_H
#define TOOLS_H

#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "matrix.h"

using namespace std;
const double y = 1.4;
const double pi = 3.14159265358979323846;

void printResults(Matrix& u, Matrix& e);

Matrix roots(double a, double b, double c);

double max(double a, double b);
double min(double a, double b);

double waveSpeed(Matrix& u, Matrix& n);

Matrix flux(Matrix& uL, Matrix& uR, Matrix& n, int R);

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
    Matrix I2E;
    Matrix B2E;
    Matrix In;
    Matrix Bn;
    Matrix Il;
    Matrix Bl;
    Matrix A;
    Matrix Br;
    Matrix Ir;
    Matrix C;

    Matrix V;
    Matrix E;
    Matrix* B;
    string* Bname;

    FVmesh(string filename);
    Matrix connectivity(Matrix& E);
    vector<int> edge2vertex(int t, int e, Matrix& E);
    Matrix boundaryConnectivity(Matrix& E, Matrix& boundary, int bgroup);
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
Matrix FV_solve(FVstate& u, FVmesh m, FVConditions c);

Matrix residual(FVstate u, FVmesh m, FVConditions c, Matrix& dt, double CFL);

void gradient(FVstate& u, FVmesh& m);

//
////Finite Element Storage
//class FEConditions 
//{
//public:
//	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//	double M, a, R, Pinf, Tt, Pt, c, p;
//	Vector2d v;
//	Vector4d uinf;
//	void set() {
//		R = 1;
//		Pinf = 1;
//		Tt = 1 + (y - 1) / 2 * M * M;
//		Pt = pow(Tt, y / (y - 1));
//		c = sqrt(y * R * Tt);
//		p = Pt / (R * Tt);
//
//
//		v(0) = M * c * cos(a);
//		v(1) = M * c * sin(a);
//
//		uinf(0) = p;
//		uinf(1) = p * v(0);
//		uinf(2) = p * v(1);
//		uinf(3) = Pt / (y - 1) + .5 * p * v.dot(v);
//		return;
//	};
//
//	FEConditions() {
//		M = 0;
//		a = 0;
//		set();
//	};
//
//	FEConditions(double Mach, double aoa) {
//		M = Mach;
//		a = aoa * pi / 180;
//		set();
//	};
//
//
//
//};
//
//
//class FEmesh
//{
//public:
//	MatrixXi I2E;
//	MatrixXi B2E;
//	MatrixXd In;
//	MatrixXd Bn;
//	MatrixXd Il;
//	MatrixXd Bl;
//	MatrixXd A;
//	MatrixXd Br;
//	MatrixXd Ir;
//	MatrixXd C;
//	Tensor nodes;
//
//	FEmesh(string filename) {
//		ifstream input;
//		int x, y, z;
//		input.open(filename);
//
//		//I2E
//		input >> x >> y;
//		I2E.resize(x, y);
//		for (int i = 0; i < x; ++i) {
//			for (int j = 0; j < y; ++j) {
//				input >> I2E(i, j);
//			}
//		}
//
//		//B2E
//		input >> x >> y;
//		B2E.resize(x, y);
//		for (int i = 0; i < x; ++i) {
//			for (int j = 0; j < y; ++j) {
//				input >> B2E(i, j);
//			}
//		}
//
//		//In
//		input >> x >> y;
//		In.resize(x, y);
//		for (int i = 0; i < x; ++i) {
//			for (int j = 0; j < y; ++j) {
//				input >> In(i, j);
//			}
//		}
//
//		//Il
//		input >> x >> y;
//		Il.resize(x, y);
//		for (int i = 0; i < x; ++i) {
//			for (int j = 0; j < y; ++j) {
//				input >> Il(i, j);
//			}
//		}
//
//
//		//Bn
//		input >> x >> y;
//		Bn.resize(x, y);
//		for (int i = 0; i < x; ++i) {
//			for (int j = 0; j < y; ++j) {
//				input >> Bn(i, j);
//			}
//		}
//
//		//Bl
//		input >> x >> y;
//		Bl.resize(x, y);
//		for (int i = 0; i < x; ++i) {
//			for (int j = 0; j < y; ++j) {
//				input >> Bl(i, j);
//			}
//		}
//
//		//A
//		input >> x >> y;
//		A.resize(x, y);
//		for (int i = 0; i < x; ++i) {
//			for (int j = 0; j < y; ++j) {
//				input >> A(i, j);
//			}
//		}
//
//		//Br
//		input >> x >> y;
//		Br.resize(x, y);
//		for (int i = 0; i < x; ++i) {
//			for (int j = 0; j < y; ++j) {
//				input >> Br(i, j);
//			}
//		}
//
//		//Ir
//		input >> x >> y;
//		Ir.resize(x, y);
//		for (int i = 0; i < x; ++i) {
//			for (int j = 0; j < y; ++j) {
//				input >> Ir(i, j);
//			}
//		}
//
//		//C
//		input >> x >> y;
//		C.resize(x, y);
//		for (int i = 0; i < x; ++i) {
//			for (int j = 0; j < y; ++j) {
//				input >> C(i, j);
//			}
//		}
//
//		//nodes
//		input >> x >> y >> z;
//		nodes.resize(x, y, z);
//		for (int i = 0; i < x; ++i) {
//			for (int j = 0; j < y; ++j) {
//				for (int k = 0; k < z; ++k) {
//					input >> nodes(i, j, k);
//				}
//			}
//		}
//
//
//		input.close();
//	};
//
//};
//
//
//class FEstate
//{
//public:
//	int p;//numerical order of accuracy
//	int q;//geometry order of accuraycy
//	MatrixXd u;
//
//	FEstate() {
//		p = 1;
//		q = 1;
//	};
//
//	FEstate(int pp, int qq, const FEmesh& m, const FEConditions& c) {
//		p = pp;
//		q = qq;
//		u.resize(m.A.rows() * (p + 1) * (p + 2) / 2, 4);
//		for (int i = 0; i < u.rows(); ++i) {
//			for (int j = 0; j < u.cols(); ++j) {
//				u(i, j) = c.uinf(j);
//			}
//		}
//	};
//
//
//};
//
//
////Finite Element Methods
//VectorXd FE_solve(FEstate& u, const FEmesh& m, const FEConditions& c);
//
//MatrixXd residual(FEstate& u, const FEmesh& m, const FEConditions& c, MatrixXd& dt, double CFL);


#endif