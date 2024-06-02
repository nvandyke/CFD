#include "Tools.h"
#include "matrix.h"
#include <assert.h>
#include <map>

//print the current state and error history to files
void printResults(Matrix& u, Matrix& e) {

    //output state
    std::ofstream output1;
    output1.open("u.txt");
    output1.precision(15);
    output1 << std::fixed;
    output1 << u << std::endl;
    output1.close();

    //output residual history
    std::ofstream output2;
    output2.open("e.txt");
    output2 << e << std::endl;
    output2.close();
}


//calculate real roots to quadratic functions
Matrix roots(double a, double b, double c) {


    double discriminant = b * b - 4 * a * c;
    Matrix ans(2, 1);
    assert(discriminant >= 0);
    double discriminant_sqrt = sqrt(discriminant);
    ans(0, 0) = (-b + discriminant_sqrt) / (2 * a);
    ans(1, 0) = (-b - discriminant_sqrt) / (2 * a);

    return ans;
}


//calculate the wave speed of the state normal to the boundary
double waveSpeed(Matrix& u, Matrix& n) {

    //velocity
    Matrix v = u.getBlock(0, 1, 1, 2) / u(0, 0);

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


//calculate the flux across the boundary uL to uR
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
    case HLLE:
    {
        //HLLE
        //HLLE flux calcs
        double smax = std::max(std::max(0.0, vL.dotProduct(vL, n) + cL), std::max(0.0, vR.dotProduct(vR, n) + cR));
        double smin = std::min(std::min(0.0, vL.dotProduct(vL, n) - cL), std::min(0.0, vR.dotProduct(vR, n) - cR));


        Matrix fHLLE(4, 1);
        fHLLE = (FuLdn + FuRdn) * .5 - (FuRdn - FuLdn) * (.5 * (smax + smin) / (smax - smin)) + (uR.transpose() - uL.transpose()) * (smax * smin / (smax - smin));


        return fHLLE;
        break;
    }
    case Rusanov:
    {
        //Rusanov
        //Rusanov flux max wave speed
        double smaxRus = std::max(sL, sR);


        Matrix fRus(4, 1);
        fRus = (FuLdn + FuRdn) * .5 - (uR.transpose() - uL.transpose()) * (.5 * smaxRus);
        return fRus;
        break;
    }
    case Roe:
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


//flux state
Matrix F(Matrix& u) {

    Matrix v = u.getBlock(0, 1, 1, 2) / u(0, 0);

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


//construct the vitrual outflow state on a subsonic outlet boundary
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


//construct the vitrual inflow state on a subsonic inlet boundary
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


//compute the flux against an inviscid wall boundary
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


//default State constructor
State::State() {
    Order = 1;
    Order_geom = 1;
    M = 0;
    a = 0;
    set();
}


//fully-defined State contructor
State::State(int O, Mesh m, double Mach, double aoa) {
    Order = O;
    mesh = m;
    M = Mach;
    a = aoa * pi / 180;
    set();

    Order_geom = 1;

    u = Matrix(mesh.A.rows(), 4);
    gradux = Matrix(mesh.A.rows(), 4);
    graduy = Matrix(mesh.A.rows(), 4);
    for (int i = 0; i < u.rows(); ++i) {
        for (int j = 0; j < u.cols(); ++j) {
            u(i, j) = uinf(0, j);
        }
    }
}


//default Mesh constructor from mesh file
Mesh::Mesh() {

}

//fully-defined Mesh constructor from mesh file
Mesh::Mesh(std::string filename) {
    std::ifstream input;
    input.open(filename);
    double double_dummy;
    int int_dummy;
    int Nn, Ne, dim;

    input >> Nn >> Ne >> dim;
    V = Matrix(Nn, dim);
    //read vertices
    for (int i = 0; i < Nn; ++i) {
        for (int j = 0; j < dim; ++j) {
            input >> double_dummy;
            V(i, j) = double_dummy;
        }
    }


    //number of boundaries
    int NB;

    input >> NB;

    //allocate space to store the info

    std::vector<Block> B;
    std::vector<std::string> Bname;

    int nBFace, nf, totalBFace = 0;
    std::string Title;

    for (int NBi = 0; NBi < NB; NBi++) {
        input >> nBFace >> nf >> Title;
        B.push_back(Block(nBFace, nf));
        Bname.push_back(Title);
        for (int j = 0; j < nBFace; ++j) {
            for (int k = 0; k < nf; ++k) {
                input >> int_dummy;
                B[NBi](j, k) = int_dummy;
            }
        }
        totalBFace += nBFace;

    }

    int Ne0 = 0;

    input >> int_dummy >> int_dummy >> Title;



    E = Block(Ne, dim + 1);
    for (int i = 0; i < Ne; ++i) {
        for (int j = 0; j < dim + 1; ++j) {
            input >> int_dummy;
            E(i, j) = int_dummy;
        }
    }

    input.close();


    //I2E
    I2E = connectivity();

    int b2erow = 0;


    B2E = Block(totalBFace, 3);
    Block B2Ei = Block(1, 3);
    int wallType;
    for (int i = 0; i < NB; i++) {
        //std::cout << Bname[i] << std::endl;
        if (Bname[i] == "Inviscid_Wall") {
            wallType = Inviscid_Wall;
        } else if (Bname[i] == "Subsonic_Outlet") {
            wallType = Subsonic_Outlet;
        } else if (Bname[i] == "Subsonic_Inlet") {
            wallType = Subsonic_Inlet;
        } else if (Bname[i] == "Supersonic_Outlet") {
            wallType = Supersonic_Outlet;
        } else if (Bname[i] == "Supersonic_Inlet") {
            wallType = Supersonic_Inlet;
        } else {
            wallType = Freestream;
        }


        B2Ei = boundaryConnectivity(B[i], wallType);
        //B2Ei = boundaryConnectivity(B[i], Freestream);
        B2E.setBlock(b2erow, 0, B2Ei);
        b2erow += B2Ei.rows();
    }


    C = Matrix(E.rows(), 2);
    In = Matrix(I2E.rows(), 2);
    Il = Matrix(I2E.rows(), 1);
    Ir = Matrix(I2E.rows(), 2);

    for (int i = 0; i < I2E.rows(); i++) {
        Matrix norm = normal(I2E(i, 0), I2E(i, 1));
        In.setBlock(i, 0, norm);
        Il(i, 0) = length(I2E(i, 0), I2E(i, 1));
        Matrix midp = midpoint(I2E(i, 0), I2E(i, 1));
        Ir.setBlock(i, 0, midp);


        Matrix temp1 = C.getBlock(I2E(i, 0), 0, 1, 2) + Ir.getBlock(i, 0, 1, 2);
        Matrix temp3 = C.getBlock(I2E(i, 2), 0, 1, 2) + Ir.getBlock(i, 0, 1, 2);

        C.setBlock(I2E(i, 0), 0, temp1);
        C.setBlock(I2E(i, 2), 0, temp3);
    }

    Bn = Matrix(B2E.rows(), 2);
    Bl = Matrix(B2E.rows(), 1);
    Br = Matrix(B2E.rows(), 2);

    for (int i = 0; i < B2E.rows(); i++) {
        Matrix norm = normal(B2E(i, 0), B2E(i, 1));
        Bn.setBlock(i, 0, norm);
        Bl(i, 0) = length(B2E(i, 0), B2E(i, 1));
        Matrix midp = midpoint(B2E(i, 0), B2E(i, 1));
        Br.setBlock(i, 0, midp);

        Matrix temp = C.getBlock(B2E(i, 0), 0, 1, 2) + Br.getBlock(i, 0, 1, 2);

        C.setBlock(B2E(i, 0), 0, temp);
    }

    C /= 3;

    A = Matrix(E.rows(), 4);
    for (int i = 0; i < E.rows(); i++) {
        A(i, 0) = 1 / Area(i);
        for (int j = 1; j < A.cols(); ++j) {
            A(i, j) = A(i, 0);
        }
    }


    std::cout << filename << " error: " << verify() << std::endl << std::endl;
    /*
    std::cout << "V" << std::endl << V << std::endl << std::endl;
    std::cout << "E" << std::endl << E << std::endl << std::endl;
    std::cout << "I2E" << std::endl << I2E << std::endl << std::endl;
    std::cout << "B2E" << std::endl << B2E << std::endl << std::endl;
    std::cout << "In" << std::endl << In << std::endl << std::endl;
    std::cout << "Il" << std::endl << Il << std::endl << std::endl;
    std::cout << "Ir" << std::endl << Ir << std::endl << std::endl;
    std::cout << "Bn" << std::endl << Bn << std::endl << std::endl;
    std::cout << "Bl" << std::endl << Bl << std::endl << std::endl;
    std::cout << "Br" << std::endl << Br << std::endl << std::endl;
    std::cout << "A" << std::endl << A << std::endl << std::endl;
    std::cout << "C" << std::endl << C << std::endl << std::endl;
    */
};


//convert triangle and edge index into the node 1/node 2 ordered pair
std::vector<int> Mesh::edge2vertex(int t, int e) {

    int n1, n2;
    std::vector<int> returnVal;

    switch (e) {
    case (0):
        if (E(t, 1) > E(t, 2)) {
            n1 = E(t, 2);
            n2 = E(t, 1);
        } else {
            n1 = E(t, 1);
            n2 = E(t, 2);
        }
        break;

    case (1):
        if (E(t, 0) > E(t, 2)) {
            n1 = E(t, 2);
            n2 = E(t, 0);
        } else {
            n1 = E(t, 0);
            n2 = E(t, 2);
        }
        break;

    case(2):
        if (E(t, 0) > E(t, 1)) {
            n1 = E(t, 1);
            n2 = E(t, 0);
        } else {
            n1 = E(t, 0);
            n2 = E(t, 1);
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


//determine the triangle/edge pairs that match up for an internal edge
Block Mesh::connectivity() {
    std::map<std::pair<int, int>, int> S;

    int n1, n2;
    int count = 0;
    std::vector<std::vector<int>> D;
    std::pair<int, int> key;
    std::pair<int, int> keyInv;

    for (int t = 0; t < E.rows(); ++t) {
        for (int e = 0; e < 3; ++e) {
            std::vector<int> nodes = edge2vertex(t, e);
            n1 = nodes[0];
            n2 = nodes[1];
            key = std::make_pair(n1, n2);
            keyInv = std::make_pair(n2, n1);
            if (S[key] > 0) {
                std::vector<int> Di;
                Di.push_back(S[key] - 1);
                Di.push_back(t);
                Di.push_back(S[keyInv] - 1);
                Di.push_back(e);
                D.push_back(Di);

                S[key] = 0;
                S[keyInv] = 0;
                count++;
            } else {
                S[key] = t + 1;
                S[keyInv] = e + 1;
            }
        }
    }

    Block C(count, 4);
    for (int i = 0; i < count; ++i) {
        C(i, 1) = D[i][2];
        C(i, 2) = D[i][1];
        C(i, 0) = D[i][0];
        C(i, 3) = D[i][3];
    }

    return C;
};


//determine the triangle/edge that match up for an boundary edge
Block Mesh::boundaryConnectivity(Block& boundary, int bgroup) {

    Block B2E(boundary.rows(), 3);
    int B2Ei = 0;

    for (int i = 0; i < boundary.rows(); i++) {

        int n1 = int(boundary(i, 0));
        int n2 = int(boundary(i, 1));
        std::vector<std::vector<int>> saved;

        for (int j = 0; j < E.rows(); j++) {
            for (int k = 0; k < 3; k++) {
                if (E(j, k) == n1) {
                    std::vector<int> here;
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
                        std::cout << "You fucked up: " << bgroup << std::endl;
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


//sets the initial conditions
void State::set() {
    Rinf = 1;
    Pinf = 1;
    if (M < 1) {
        Tt = 1 + (y - 1) / 2 * M * M;
        Pt = pow(Tt, y / (y - 1));
        c = sqrt(y * Rinf * Tt);
        p = Pt / (Rinf * Tt);

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


//convert triangle and edge index into the XY points
Matrix Mesh::pointsFromTE(int t, int e) {
    int n1, n2;
    switch (e) {
    case 0:
        n1 = E(t, 1) - 1;
        n2 = E(t, 2) - 1;
        break;
    case 1:
        n1 = E(t, 2) - 1;
        n2 = E(t, 0) - 1;
        break;
    case 2:
        n1 = E(t, 0) - 1;
        n2 = E(t, 1) - 1;
        break;
    default:
        std::cout << "poop" << std::endl;
        assert(false);
        break;
    }
    Matrix pts(3, 2);
    pts(0, 0) = V(n1, 0);
    pts(0, 1) = V(n1, 1);
    pts(1, 0) = V(n2, 0);
    pts(1, 1) = V(n2, 1);
    pts(2, 0) = V(E(t, e) - 1, 0);
    pts(2, 1) = V(E(t, e) - 1, 1);

    return pts;

}


//get the normal vector from triangle and edge index
Matrix Mesh::normal(int t, int e) {
    Matrix pts = pointsFromTE(t, e);

    Matrix normal(1, 2);
    normal(0, 0) = pts(1, 1) - pts(0, 1);
    normal(0, 1) = pts(0, 0) - pts(1, 0);
    double l = length(t, e);
    normal /= l;
    return normal;
}


//get the length of the edge from the triangle and edge index
double Mesh::length(int t, int e) {
    Matrix pts = pointsFromTE(t, e);
    double length = sqrt((pts(0, 0) - pts(1, 0)) * (pts(0, 0) - pts(1, 0)) + (pts(0, 1) - pts(1, 1)) * (pts(0, 1) - pts(1, 1)));
    return length;
}


//get the midpoint of the edge from the triangle and edge index
Matrix Mesh::midpoint(int t, int e) {
    Matrix pts = pointsFromTE(t, e);

    Matrix midpoint(1, 2);
    midpoint(0, 0) = (pts(0, 0) + pts(1, 0)) / 2.0;
    midpoint(0, 1) = (pts(0, 1) + pts(1, 1)) / 2.0;

    return midpoint;

}


//calculate the area of the triangle
double Mesh::Area(int t) {
    Matrix pts = pointsFromTE(t, 0);

    double A = abs(.5 * (pts(0, 0) * (pts(1, 1) - pts(2, 1)) + pts(1, 0) * (pts(2, 1) - pts(0, 1)) + pts(2, 0) * (pts(0, 1) - pts(1, 1))));
    assert(A > 1e-10);
    return A;
}


//verify the mesh is imported correctly, there should be zero residual
double Mesh::verify() {
    Matrix sum(E.rows(), 2);
    Matrix tot(E.rows(), 1);
    Matrix curSum(1, 2);
    Matrix curSum1(1, 2);
    Matrix curSum2(1, 2);

    for (int i = 0; i < I2E.rows(); ++i) {
        curSum1 = sum.getBlock(I2E(i, 0), 0, 1, 2) + In.getBlock(i, 0, 1, 2) * Il(i, 0);
        curSum2 = sum.getBlock(I2E(i, 2), 0, 1, 2) - In.getBlock(i, 0, 1, 2) * Il(i, 0);

        sum.setBlock(I2E(i, 0), 0, curSum1);
        sum.setBlock(I2E(i, 2), 0, curSum2);
    }
    for (int i = 0; i < B2E.rows(); ++i) {
        curSum = sum.getBlock(B2E(i, 0), 0, 1, 2) + Bn.getBlock(i, 0, 1, 2) * Bl(i, 0);

        sum.setBlock(B2E(i, 0), 0, curSum);
    }
    for (int i = 0; i < tot.rows(); ++i) {
        tot(i, 0) = sqrt(sum(i, 0) * sum(i, 0) + sum(i, 1) * sum(i, 1));
    }
    return tot.max();

}
