#include "Tools.h"
#include "matrix.h"
#include <assert.h>

void printResults(Matrix& u, Matrix& e) {

    //output state
    ofstream output1;
    output1.open("u.txt");
    output1.precision(15);
    output1 << std::fixed;
    output1 << u << endl;
    output1.close();

    //output residual history
    ofstream output2;
    output2.open("e.txt");
    output2 << e << endl;
    output2.close();
}


Matrix roots(double a, double b, double c) {
    //calculate real roots to quadratic functions

    double discriminant = b * b - 4 * a * c;
    Matrix ans(2, 1);
    assert(discriminant >= 0);
    double discriminant_sqrt = sqrt(discriminant);
    ans(0, 0) = (-b + discriminant_sqrt) / (2 * a);
    ans(1, 0) = (-b - discriminant_sqrt) / (2 * a);

    return ans;
}


double max(double a, double b) {
    if (a > b) {
        return a;
    } else {
        return b;
    }
}

double min(double a, double b) {
    if (a < b) {
        return a;
    } else {
        return b;
    }
}

double waveSpeed(Matrix& u, Matrix& n) {

    //velocity
    Matrix v = u.getBlock(0, 1, 1, 2);
    v(0, 0) = v(0, 0) / u(0, 0);
    v(0, 1) = v(0, 1) / u(0, 0);

    //pressure
    double P = (y - 1) * (u(0, 3) - .5 * u(0, 0) * v.dotProduct(v, v));

    //speed of sound
    double c = sqrt(y * P / u(0, 0));

    //normal velocity
    double vn = v.dotProduct(v, n);

    //wave speed
    double s = fabs(vn) + c;
    return s;
}


Matrix flux(Matrix& uL, Matrix& uR, Matrix& n, int R) {

    //velocities
    Matrix vL = uL.getBlock(0, 1, 1, 2);
    Matrix vR = uR.getBlock(0, 1, 1, 2);
    Matrix vdiff = vR - vL;

    //Left state
    double rhoL = uL(0, 0);
    double rhoEL = uL(0, 3);
    vL(0, 0) = vL(0, 0) / rhoL;
    vL(0, 1) = vL(0, 1) / rhoL;
    double halfvLdotvL = .5 * vL.dotProduct(vL, vL);
    double PL = (y - 1) * (rhoEL - rhoL * halfvLdotvL);
    double HL = (rhoEL + PL) / rhoL;
    double cL = (y - 1) * (HL - halfvLdotvL);
    double sL = abs(vL.dotProduct(vL, n)) + cL;

    //Right state
    double rhoR = uR(0, 0);
    double rhoER = uR(0, 3);
    vR(0, 0) = vR(0, 0) / rhoR;
    vR(0, 1) = vR(0, 1) / rhoR;
    double halfvRdotvR = .5 * vR.dotProduct(vR, vR);
    double PR = (y - 1) * (rhoER - rhoR * halfvRdotvR);
    double HR = (rhoER + PR) / rhoR;
    double cR = (y - 1) * (HR - halfvRdotvR);
    double sR = abs(vR.dotProduct(vR, n)) + cR;

    //numerical flux calculation
    Matrix nT = n.transpose();
    Matrix FuLdn = F(uL) * nT;
    Matrix FuRdn = F(uR) * nT;

    switch (R) {
    case 1:
    {
        //HLLE
        //HLLE flux calcs
        double smax = max(max(0.0, vL.dotProduct(vL, n) + cL), max(0.0, vR.dotProduct(vR, n) + cR));
        double smin = min(min(0.0, vL.dotProduct(vL, n) - cL), min(0.0, vR.dotProduct(vR, n) - cR));


        Matrix fHLLE(4, 1);
        fHLLE = (FuLdn + FuRdn) * .5 - (FuRdn - FuLdn) * (.5 * (smax + smin) / (smax - smin)) + (uR.transpose() - uL.transpose()) * (smax * smin / (smax - smin));


        return fHLLE;
        break;
    }
    case 2:
    {
        //Rusanov
        //Rusanov flux max wave speed
        double smaxRus = max(sL, sR);


        Matrix fRus(4, 1);
        fRus = (FuLdn + FuRdn) * .5 - (uR.transpose() - uL.transpose()) * (.5 * smaxRus);
        return fRus;
        break;
    }
    case 0:
    default:
    {
        //Roe

        //Roe-averaged state
        //Matrix v(1, 2);
        double sqrt_rhoL = sqrt(rhoL);
        double sqrt_rhoR = sqrt(rhoR);
        double denom = 1 / (sqrt_rhoL + sqrt_rhoR);


        Matrix v = (vL * sqrt_rhoL + vR * sqrt_rhoR) * denom;
        double H = (sqrt_rhoL * HL + sqrt_rhoR * HR) * denom;
        double u = v.dotProduct(v, n);
        double q2 = v.dotProduct(v, v);
        double c = (y - 1) * (H - .5 * q2);

        //eigenvalues
        double lambda[3] = { abs(u + c), abs(u - c), u };

        double e = 0.1 * c;
        for (int i = 0; i < 3; ++i) {
            if (lambda[i] < e)
                lambda[i] = abs((e * e + lambda[i] * lambda[i]) / (2 * e));
        }

        //subcalcs for A matrix
        double s1 = .5 * (lambda[0] + lambda[1]) - lambda[2];
        double s2 = .5 * (lambda[0] - lambda[1]);

        double G1 = ((y - 1) * (q2 / 2 * (rhoR - rhoL) - v.dotProduct(v, vdiff) + (rhoER - rhoEL))) / c;
        double G2 = -u * (rhoR - rhoL) + n.dotProduct(n, vdiff);

        double C1 = (G1 * s1 + G2 * s2) / c;
        double C2 = G1 * s2 + s1 * G2;

        //upwinding
        Matrix A(4, 1);
        A(0, 0) = lambda[2] * (rhoR - rhoL) + C1;
        A(1, 0) = lambda[2] * (vdiff(0, 0)) + C1 * v(0, 0) + C2 * n(0, 0);
        A(2, 0) = lambda[2] * (vdiff(0, 1)) + C1 * v(0, 1) + C2 * n(0, 1);
        A(3, 0) = lambda[2] * (rhoER - rhoEL) + C1 * H + C2 * u;


        Matrix fRoe = (FuLdn + FuRdn - A) * .5;
        return fRoe;
        break;
    }
    }

}


Matrix F(Matrix& u) {

    Matrix v = u.getBlock(0, 1, 1, 2);
    v(0, 0) = v(0, 0) / u(0, 0);
    v(0, 1) = v(0, 1) / u(0, 0);

    double kinematicenergy = (y - 1) * (u(0, 3) - (0.5 * u(0, 0) * v.dotProduct(v, v)));

    Matrix f(4, 2);
    f(0, 0) = u(0, 1);
    f(0, 1) = u(0, 2);

    f(1, 0) = u(0, 1) * v(0, 0) + kinematicenergy;
    f(1, 1) = u(0, 1) * v(0, 1);

    f(2, 0) = u(0, 2) * v(0, 0);
    f(2, 1) = u(0, 2) * v(0, 1) + kinematicenergy;

    f(3, 0) = u(0, 3) * v(0, 0) + v(0, 0) * kinematicenergy;
    f(3, 1) = u(0, 3) * v(0, 1) + v(0, 1) * kinematicenergy;


    return f;
}


Matrix Outflow(Matrix& u, Matrix& n, double Pb) {


    //internal state quantities

    //velocity
    Matrix v = u.getBlock(0, 1, 1, 2);
    v(0, 0) = v(0, 0) / u(0, 0);
    v(0, 1) = v(0, 1) / u(0, 0);

    //pressure
    double P = (y - 1) * (u(0, 3) - .5 * u(0, 0) * v.dotProduct(v, v));

    //speed of sound
    double c = sqrt(y * P / u(0, 0));

    //entropy
    double S = P / pow(u(0, 0), y);

    //boundary state quantities

    //density
    double pb = pow(Pb / S, 1 / y);

    //speed of sound
    double cb = sqrt(y * Pb / pb);

    //normal component of velocity
    double ubn = v.dotProduct(v, n) + 2 * (c - cb) / (y - 1);

    //velocity
    Matrix vb = v - n * v.dotProduct(v, n) + n * ubn;

    //energy term
    double pEb = Pb / (y - 1) + .5 * pb * vb.dotProduct(vb, vb);


    //state construction
    Matrix ub(1, 4);
    ub(0, 0) = pb;
    ub(0, 1) = pb * vb(0, 0);
    ub(0, 2) = pb * vb(0, 1);
    ub(0, 3) = pEb;
    return ub;
}


Matrix Inflow(Matrix& u, Matrix& n, double Tt, double Pt, double a, double R) {


    //internal state quantities

    //velocity
    Matrix v = u.getBlock(0, 1, 1, 2);
    v(0, 0) = v(0, 0) / u(0, 0);
    v(0, 1) = v(0, 1) / u(0, 0);

    //pressure
    double P = (y - 1) * (u(0, 3) - .5 * u(0, 0) * v.dotProduct(v, v));

    //speed of sound
    double c = sqrt(y * P / u(0, 0));

    //riemann invariant J+
    double J = v.dotProduct(v, n) + 2 * c / (y - 1);


    //quantities for computing Mach #
    Matrix nin(1, 2);
    nin(0, 0) = cos(a);
    nin(0, 1) = sin(a);
    double dn = nin.dotProduct(nin, n);

    Matrix Mb = roots(y * R * Tt * dn * dn - (y - 1) / 2 * J * J, 4 * y * R * Tt * dn / (y - 1), 4 * y * R * Tt / (y - 1) / (y - 1) - J * J);
    double MB;

    //since it's quadratic, we need to grab appropriate number
    if (Mb(0, 0) > 0 && Mb(1, 0) > 0) {
        MB = Mb.max();
    } else if (Mb(0, 0) > 0) {
        MB = Mb(0, 0);
    } else if (Mb(1, 0) > 0) {
        MB = Mb(1, 0);
    } else {
        //catch errors
        std::cout << "you fucked up, ";
        (Mb.transpose()).print();
        u.print();
        n.print();
        MB = 0.5;
    }

    //define boundary state
    double Tb = Tt / (1 + .5 * (y - 1) * MB * MB);
    double Pb = Pt * pow(Tb / Tt, y / (y - 1));
    double pb = Pb / (R * Tb);
    double cb = sqrt(y * Pb / pb);
    Matrix vb = nin * (MB * cb);
    double pEb = Pb / (y - 1) + .5 * pb * vb.dotProduct(vb, vb);

    //output
    Matrix ub(1, 4);
    ub(0, 0) = pb;
    ub(0, 1) = pb * vb(0, 0);
    ub(0, 2) = pb * vb(0, 1);
    ub(0, 3) = pEb;
    return ub;
}


Matrix wallFlux(Matrix& u, Matrix& n) {

    //velocity
    Matrix v = u.getBlock(0, 1, 1, 2);
    v(0, 0) = v(0, 0) / u(0, 0);
    v(0, 1) = v(0, 1) / u(0, 0);

    //velocity at boundary
    Matrix vb = v - n * v.dotProduct(v, n);

    //pressure at boundary
    double Pb = (y - 1) * (u(0, 3) - .5 * u(0, 0) * v.dotProduct(vb, vb));

    //flux
    Matrix f(1, 4);
    f(0, 0) = 0;
    f(0, 1) = Pb * n(0, 0);
    f(0, 2) = Pb * n(0, 1);
    f(0, 3) = 0;
    return f;
}



FVstate::FVstate() {
    Order = 1;
}
FVstate::FVstate(int O, FVmesh m, FVConditions c) {
    Order = O;
    u = Matrix(m.A.rows(), 4);
    gradux = Matrix(m.A.rows(), 4);
    graduy = Matrix(m.A.rows(), 4);
    for (int i = 0; i < u.rows(); ++i) {
        for (int j = 0; j < u.cols(); ++j) {
            u(i, j) = c.uinf(0, j);
        }
    }
};


FVmesh::FVmesh(string filename) {
    ifstream input;
    input.open(filename);
    double dummy;
    int Nn, Ne, dim;

    input >> Nn >> Ne >> dim;
    V = Matrix(Nn, dim);
    //read vertices
    for (int i = 0; i < Nn; ++i) {
        for (int j = 0; j < dim; ++j) {
            input >> dummy;
            V(i, j) = dummy;
        }
    }


    //number of boundaries
    int NB;

    input >> NB;

    //allocate space to store the info, 
    B = new Matrix[NB];
    Bname = new string[NB];

    int nBFace, nf, totalBFace = 0;
    string Title;

    for (int NBi = 0; NBi < NB; NBi++) {
        input >> nBFace >> nf >> Title;
        B[NBi] = Matrix(nBFace, nf);
        Bname[NBi] = Title;
        for (int j = 0; j < nBFace; ++j) {
            for (int k = 0; k < nf; ++k) {
                input >> dummy;
                B[NBi](j, k) = dummy;
            }
        }
        totalBFace += nBFace;

    }

    int Ne0 = 0;

    input >> dummy >> dummy >> Title;



    E = Matrix(Ne, dim + 1);
    for (int i = 0; i < Ne; ++i) {
        for (int j = 0; j < dim + 1; ++j) {
            input >> dummy;
            E(i, j) = dummy;
        }
    }

    input.close();


    //I2E
    I2E = connectivity(E);

    int b2erow = 0;


    B2E = Matrix(totalBFace, 3);
    for (int i = 0; i < NB; i++) {
        //swap for freestream test
        Matrix B2Ei = boundaryConnectivity(E, B[i], i + 1);
        //Matrix B2Ei = boundaryConnectivity(E, B[i], 7);
        B2E.setBlock(b2erow, 0, B2Ei);
        b2erow += B2Ei.rows();
    }


    C = Matrix(E.rows(), 2);
    In = Matrix(I2E.rows(), 2);
    Il = Matrix(I2E.rows(), 1);
    Ir = Matrix(I2E.rows(), 2);

    for (int i = 0; i < I2E.rows(); i++) {
        Matrix norm = normal(int(I2E(i, 0)), int(I2E(i, 1)));
        In.setBlock(i, 0, norm);
        Il(i, 0) = length(int(I2E(i, 0)), int(I2E(i, 1)));
        Matrix midp = midpoint(int(I2E(i, 0)), int(I2E(i, 1)));
        Ir.setBlock(i, 0, midp);


        Matrix temp1 = C.getBlock(int(I2E(i, 0)), 0, 1, 2) + Ir.getBlock(i, 0, 1, 2);
        Matrix temp3 = C.getBlock(int(I2E(i, 2)), 0, 1, 2) + Ir.getBlock(i, 0, 1, 2);

        C.setBlock(int(I2E(i, 0)), 0, temp1);
        C.setBlock(int(I2E(i, 2)), 0, temp3);
    }

    Bn = Matrix(B2E.rows(), 2);
    Bl = Matrix(B2E.rows(), 1);
    Br = Matrix(B2E.rows(), 2);

    for (int i = 0; i < B2E.rows(); i++) {
        Matrix norm = normal(int(B2E(i, 0)), int(B2E(i, 1)));
        Bn.setBlock(i, 0, norm);
        Bl(i, 0) = length(int(B2E(i, 0)), int(B2E(i, 1)));
        Matrix midp = midpoint(int(B2E(i, 0)), int(B2E(i, 1)));
        Br.setBlock(i, 0, midp);

        Matrix temp = C.getBlock(int(B2E(i, 0)), 0, 1, 2) + Br.getBlock(i, 0, 1, 2);

        C.setBlock(int(B2E(i, 0)), 0, temp);
    }

    C /= 3;

    A = Matrix(E.rows(), 4);
    for (int i = 0; i < E.rows(); i++) {
        A(i, 0) = 1 / Area(i);
        for (int j = 1; j < A.cols(); ++j) {
            A(i, j) = A(i, 0);
        }
    }


    cout << filename << " error: " << verify() << endl << endl;
    /*
    cout << "V" << endl << V << endl << endl;
    cout << "E" << endl << E << endl << endl;
    cout << "I2E" << endl << I2E << endl << endl;
    cout << "B2E" << endl << B2E << endl << endl;
    cout << "In" << endl << In << endl << endl;
    cout << "Il" << endl << Il << endl << endl;
    cout << "Ir" << endl << Ir << endl << endl;
    cout << "Bn" << endl << Bn << endl << endl;
    cout << "Bl" << endl << Bl << endl << endl;
    cout << "Br" << endl << Br << endl << endl;
    cout << "A" << endl << A << endl << endl;
    cout << "C" << endl << C << endl << endl;
    */
};


vector<int> FVmesh::edge2vertex(int t, int e, Matrix& E) {

    int n1, n2;
    vector<int> returnVal;

    switch (e) {
    case (0):
        if (E(t, 1) > E(t, 2)) {
            n1 = int(E(t, 2));
            n2 = int(E(t, 1));
        } else {
            n1 = int(E(t, 1));
            n2 = int(E(t, 2));
        }
        break;

    case (1):
        if (E(t, 0) > E(t, 2)) {
            n1 = int(E(t, 2));
            n2 = int(E(t, 0));
        } else {
            n1 = int(E(t, 0));
            n2 = int(E(t, 2));
        }
        break;

    case(2):
        if (E(t, 0) > E(t, 1)) {
            n1 = int(E(t, 1));
            n2 = int(E(t, 0));
        } else {
            n1 = int(E(t, 0));
            n2 = int(E(t, 1));
        }
        break;

    default:
        n1 = 0;
        n2 = 0;
        std::cout << "ERROR in edge2vertex\n";
        break;

    }
    returnVal.push_back(n1 - 1);
    returnVal.push_back(n2 - 1);
    return returnVal;
};


Matrix FVmesh::connectivity(Matrix& E) {
    Matrix S(int(E.max()), int(E.max()));

    int n1, n2;
    int count = 0;
    vector<vector<int>> D;

    for (int t = 0; t < E.rows(); ++t) {
        for (int e = 0; e < 3; ++e) {
            vector<int> nodes = edge2vertex(t, e, E);
            n1 = nodes[0];
            n2 = nodes[1];
            if (S(n1, n2) > 0) {
                vector<int> Di;
                Di.push_back(int(S(n1, n2)) - 1);
                Di.push_back(t);
                Di.push_back(int(S(n2, n1)) - 1);
                Di.push_back(e);
                D.push_back(Di);

                S(n1, n2) = 0;
                S(n2, n1) = 0;
                count++;
            } else {
                S(n1, n2) = t + 1;
                S(n2, n1) = e + 1;
            }
        }
    }

    Matrix C(count, 4);
    for (int i = 0; i < count; ++i) {
        C(i, 1) = D[i][2];
        C(i, 2) = D[i][1];
        C(i, 0) = D[i][0];
        C(i, 3) = D[i][3];
    }

    return C;
};

Matrix FVmesh::boundaryConnectivity(Matrix& E, Matrix& boundary, int bgroup) {

    Matrix B2E(boundary.rows(), 3);
    int B2Ei = 0;

    for (int i = 0; i < boundary.rows(); i++) {

        int n1 = int(boundary(i, 0));
        int n2 = int(boundary(i, 1));
        vector<vector<int>> saved;

        for (int j = 0; j < E.rows(); j++) {
            for (int k = 0; k < 3; k++) {
                if (E(j, k) == n1) {
                    vector<int> here;
                    here.push_back(j);
                    here.push_back(k);
                    saved.push_back(here);

                    break;
                }
            }
        }

        for (int j = 0; j < saved.size(); j++) {
            for (int k = 0; k < 3; k++) {

                if (E(saved[j][0], k) == n2) {
                    int face = -1;
                    if (saved[j][1] == 0 and k == 1 or saved[j][1] == 1 and k == 0) {
                        face = 2;
                    } else if (saved[j][1] == 1 and k == 2 or saved[j][1] == 2 and k == 1) {
                        face = 0;
                    } else if (saved[j][1] == 2 and k == 0 or saved[j][1] == 0 and k == 2) {
                        face = 1;
                    } else {
                        cout << "You fucked up: " << bgroup << endl;
                    }

                    B2E(B2Ei, 0) = saved[j][0];
                    B2E(B2Ei, 1) = face;
                    B2E(B2Ei, 2) = bgroup;
                    B2Ei++;
                    break;

                }
            }
        }

    }

    return B2E;
}




void FVConditions::set() {
    R = 1;
    Pinf = 1;
    if (M < 1) {
        Tt = 1 + (y - 1) / 2 * M * M;
        Pt = pow(Tt, y / (y - 1));
        c = sqrt(y * R * Tt);
        p = Pt / (R * Tt);

        v = Matrix(1, 2);
        v(0, 0) = M * c * cos(a);
        v(0, 1) = M * c * sin(a);

        uinf = Matrix(1, 4);
        uinf(0, 0) = p;
        uinf(0, 1) = p * v(0, 0);
        uinf(0, 2) = p * v(0, 1);
        uinf(0, 3) = Pt / (y - 1) + .5 * p * v.dotProduct(v, v);
    } else {
        uinf = Matrix(1, 4);
        uinf(0, 0) = 1;
        uinf(0, 1) = M * cos(a);
        uinf(0, 2) = M * sin(a);
        uinf(0, 3) = 1 / ((y - 1) * y) + .5 * M * M;
    }
    return;
};

FVConditions::FVConditions() {
    M = 0;
    a = 0;

    set();
};

FVConditions::FVConditions(double Mach, double aoa) {
    M = Mach;
    a = aoa * pi / 180;
    set();
};

Matrix FVmesh::pointsFromTE(int t, int e) {
    int n1, n2;
    switch (e) {
    case 0:
        n1 = int(E(t, 1)) - 1;
        n2 = int(E(t, 2)) - 1;
        break;
    case 1:
        n1 = int(E(t, 2)) - 1;
        n2 = int(E(t, 0)) - 1;
        break;
    case 2:
        n1 = int(E(t, 0)) - 1;
        n2 = int(E(t, 1)) - 1;
        break;
    default:
        cout << "poop" << endl;
        assert(false);
        break;
    }
    Matrix pts(3, 2);
    pts(0, 0) = V(n1, 0);
    pts(0, 1) = V(n1, 1);
    pts(1, 0) = V(n2, 0);
    pts(1, 1) = V(n2, 1);
    pts(2, 0) = V(int(E(t, e)) - 1, 0);
    pts(2, 1) = V(int(E(t, e)) - 1, 1);

    return pts;

}




Matrix FVmesh::normal(int t, int e) {
    Matrix pts = pointsFromTE(t, e);

    Matrix normal(1, 2);
    normal(0, 0) = pts(1, 1) - pts(0, 1);
    normal(0, 1) = pts(0, 0) - pts(1, 0);
    double l = length(t, e);
    normal /= l;
    return normal;
}


double FVmesh::length(int t, int e) {
    Matrix pts = pointsFromTE(t, e);
    double length = sqrt((pts(0, 0) - pts(1, 0)) * (pts(0, 0) - pts(1, 0)) + (pts(0, 1) - pts(1, 1)) * (pts(0, 1) - pts(1, 1)));
    return length;
}

Matrix FVmesh::midpoint(int t, int e) {
    Matrix pts = pointsFromTE(t, e);

    Matrix midpoint(1, 2);
    midpoint(0, 0) = (pts(0, 0) + pts(1, 0)) / 2.0;
    midpoint(0, 1) = (pts(0, 1) + pts(1, 1)) / 2.0;

    return midpoint;

}

double FVmesh::Area(int t) {
    Matrix pts = pointsFromTE(t, 0);

    double A = abs(.5 * (pts(0, 0) * (pts(1, 1) - pts(2, 1)) + pts(1, 0) * (pts(2, 1) - pts(0, 1)) + pts(2, 0) * (pts(0, 1) - pts(1, 1))));
    if (A < 1e-10) {
        cout << pts << endl;
        assert(0);
    }
    return A;
}


double FVmesh::verify() {
    Matrix sum(E.rows(), 2);
    Matrix tot(E.rows(), 1);

    for (int i = 0; i < I2E.rows(); i++) {
        Matrix curSum1 = sum.getBlock(int(I2E(i, 0)), 0, 1, 2);
        curSum1 += In.getBlock(i, 0, 1, 2) * Il(i, 0);

        Matrix curSum2 = sum.getBlock(int(I2E(i, 2)), 0, 1, 2);
        curSum2 -= In.getBlock(i, 0, 1, 2) * Il(i, 0);

        sum.setBlock(int(I2E(i, 0)), 0, curSum1);
        sum.setBlock(int(I2E(i, 2)), 0, curSum2);
    }
    for (int i = 0; i < B2E.rows(); i++) {
        Matrix curSum = sum.getBlock(int(B2E(i, 0)), 0, 1, 2);
        curSum += Bn.getBlock(i, 0, 1, 2) * Bl(i, 0);

        sum.setBlock(int(B2E(i, 0)), 0, curSum);
    }
    for (int i = 0; i < tot.rows(); i++) {
        tot(i, 0) = sqrt(sum(i, 0) * sum(i, 0) + sum(i, 1) * sum(i, 1));
    }
    return tot.max();

}