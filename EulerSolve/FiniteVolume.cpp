#include "Tools.h"
#include <omp.h>

const double CFL = 0.75;


extern "C" {
    void FVrun(string, string, int, double, double);
    void cudaStart(int);
    void cudaEnd();
}


void FVrun(string ic, string mesh, int order, double Mach, double angleOfAttack) {
    FVmesh m(mesh);

    FVConditions c(Mach, angleOfAttack);
    FVstate u(order, m, c);

    if (ic != "") {
        ifstream inputfile;
        inputfile.open(ic);
        inputfile >> u.u;
        inputfile.close();
    }

    Matrix e = FV_solve(u, m, c);
    printResults(u.u, e);
    return;
}


Matrix FV_solve(FVstate& u, FVmesh& m, FVConditions c) {
    //constants
    double tol = 1e-7;
    double max = 1 / tol;
    int MaxIter = 100000;

    //initialize
    cudaStart(u.u.size());
    Matrix e(MaxIter + 1, 1);
    Matrix dt(u.u.rows(), u.u.cols());
    Matrix R(u.u.rows(), u.u.cols());
    int iter = 0;

    //iterate
    if (u.Order == 1) {

        for (iter = 0; iter < MaxIter; ++iter) {
            R = residual(u, m, c, dt);
            u.u = u.u - dt % R;
            //u.u = u.u - dt.multInPlace(R);
            e(iter, 0) = R.max();
            cout << e(iter, 0) << endl;
            if (e(iter, 0) < tol)
                break;
            if (e(iter, 0) > max || isnan(e(iter, 0)))
                break;
            if (iter % 1000 == 0)
                printResults(u.u, e);
        }
    } else if (u.Order == 2) {

        FVstate ufe(2, m, c);
        Matrix dtdummy(u.u.rows(), u.u.cols());

        for (iter = 0; iter < MaxIter; ++iter) {
            R = residual(u, m, c, dt);
            ufe.u = u.u - dt % R;
            R = residual(ufe, m, c, dtdummy);
            u.u = (u.u + ufe.u - dt % R) * .5;
            e(iter, 0) = R.max();
            cout << e(iter, 0) << endl;
            if (e(iter, 0) < tol)
                break;
            if (e(iter, 0) > max || isnan(e(iter, 0)))
                break;
            if (iter % 1000 == 0)
                printResults(u.u, e);
        }
    } else if (u.Order < 5) {

        FVstate ufe(2, m, c);
        Matrix dtdummy(u.u.rows(), u.u.cols());
        Matrix R2(u.u.rows(), u.u.cols());
        Matrix R3(u.u.rows(), u.u.cols());
        Matrix R4(u.u.rows(), u.u.cols());
        Matrix bigR(u.u.rows(), u.u.cols());

        for (iter = 0; iter < MaxIter; ++iter) {
            R = residual(u, m, c, dt);
            ufe.u = u.u - .5 * (dt % R);
            R2 = residual(ufe, m, c, dtdummy);
            ufe.u = u.u - .5 * (dt % R2);
            R3 = residual(ufe, m, c, dtdummy);
            ufe.u = u.u - dt % R3;
            R4 = residual(ufe, m, c, dtdummy);
            bigR = R + 2 * (R2 + R3) + R4;
            u.u = u.u - (dt % bigR) / 6;
            e(iter, 0) = R.max();
            cout << e(iter, 0) << endl;
            if (e(iter, 0) < tol)
                break;
            if (e(iter, 0) > max || isnan(e(iter, 0)))
                break;
            if (iter % 1000 == 0)
                printResults(u.u, e);
        }
    }
    cudaEnd();
    cout << "Stopped at iteration " << iter << endl;
    Matrix retVal(int(iter) + 1, 1);
    retVal = e.getBlock(0, 0, iter + 1, 1);
    return retVal;
}

void interioredges(FVstate& u, FVmesh& m, double* R, double* sl) {


    Matrix gradux_l = Matrix(1, 4);
    Matrix graduy_l = Matrix(1, 4);
    Matrix gradux_r = Matrix(1, 4);
    Matrix graduy_r = Matrix(1, 4);
    Matrix uL = Matrix(1, 4);
    Matrix uR = Matrix(1, 4);
    Matrix R_l = Matrix(1, 4);
    Matrix R_r = Matrix(1, 4);
    Matrix f = Matrix(4, 1);
    Matrix fTl = Matrix(1, 4);
    Matrix n = Matrix(1, 2);
    int i, j, k;
    double l, ws;



#pragma omp parallel
    {
#pragma omp for private(j,k,uL,uR,gradux_l,graduy_l,gradux_r,graduy_r,n,l,f,ws,fTl)
        for (i = 0; i < m.I2E.rows(); ++i) {

            j = m.I2E(i, 0);
            k = m.I2E(i, 2);
            uL = u.u.getBlock(j, 0, 1, 4);
            uR = u.u.getBlock(k, 0, 1, 4);

            if (u.Order == 2) {
                gradux_l = u.gradux.getBlock(j, 0, 1, 4);
                graduy_l = u.graduy.getBlock(j, 0, 1, 4);
                gradux_r = u.gradux.getBlock(k, 0, 1, 4);
                graduy_r = u.graduy.getBlock(k, 0, 1, 4);

                uL += (gradux_l * (m.Ir(i, 0) - m.C(j, 0)) + graduy_l * (m.Ir(i, 1) - m.C(j, 1)));
                uR += (gradux_r * (m.Ir(i, 0) - m.C(k, 0)) + graduy_r * (m.Ir(i, 1) - m.C(k, 1)));

            }
            n = m.In.getBlock(i, 0, 1, 2);
            l = m.Il(i, 0);
            f = flux(uL, uR, n, Roe);

            ws = max(waveSpeed(uL, n), waveSpeed(uR, n)) * l;

            fTl = f.transpose() * l;
#pragma omp critical 
            {
                R[u.u.cols() * j + 0] += fTl(0, 0);
                R[u.u.cols() * j + 1] += fTl(0, 1);
                R[u.u.cols() * j + 2] += fTl(0, 2);
                R[u.u.cols() * j + 3] += fTl(0, 3);

                R[u.u.cols() * k + 0] -= fTl(0, 0);
                R[u.u.cols() * k + 1] -= fTl(0, 1);
                R[u.u.cols() * k + 2] -= fTl(0, 2);
                R[u.u.cols() * k + 3] -= fTl(0, 3);

                sl[j] += ws;
                sl[k] += ws;
            }
        }
    }
}

void boundaryedges(FVstate& u, FVmesh& m, FVConditions c, double* R, double* sl) {


    Matrix gradux_i = Matrix(1, 4);
    Matrix graduy_i = Matrix(1, 4);
    Matrix ui = Matrix(1, 4);
    Matrix uout = Matrix(1, 4);
    Matrix block = Matrix(1, 4);
    Matrix f = Matrix(4, 1);
    Matrix n = Matrix(1, 2);
    Matrix uin = Matrix(1, 4);
    int i, j;
    double l, ws;

#pragma omp parallel
    {
#pragma omp for private(j,ui,gradux_i,graduy_i,n,l,uout,uin,block,ws)
        for (i = 0; i < m.B2E.rows(); ++i) {
            j = m.B2E(i, 0);
            ui = u.u.getBlock(j, 0, 1, 4);
            if (u.Order == 2) {
                gradux_i = u.gradux.getBlock(j, 0, 1, 4);
                graduy_i = u.graduy.getBlock(j, 0, 1, 4);
                ui += (gradux_i * (m.Br(i, 0) - m.C(j, 0)) + graduy_i * (m.Br(i, 1) - m.C(j, 1)));
            }
            n = m.Bn.getBlock(i, 0, 1, 2);
            l = m.Bl(i, 0);

            switch ((int)m.B2E(i, 2)) {
            case Inviscid_Wall:
                //Inviscid Wall
                block = wallFlux(ui, n) * l;
                ws = fabs(waveSpeed(ui, n)) * l;
                break;

            case Subsonic_Outlet:
                //Subsonic Outlet
                uout = Outflow(ui, n, c.Pinf);
                block = (F(uout) * (n.transpose())).transpose() * l;
                ws = fabs(waveSpeed(uout, n)) * l;
                break;

            case Subsonic_Inlet:
                //Subsonic Inlet
                uin = Inflow(ui, n, c.Tt, c.Pt, c.a, c.R);
                block = (F(uin) * (n.transpose())).transpose() * l;
                ws = fabs(waveSpeed(uin, n)) * l;
                break;

            case Supersonic_Inlet:
                //Supersonic Inlet
                block = flux(ui, c.uinf, n, Roe).transpose() * l;
                ws = fabs(max(waveSpeed(ui, n), waveSpeed(c.uinf, n))) * l;
                break;

            case Supersonic_Outlet:
                //Supersonic Outlet
                block = (F(ui) * (n.transpose())).transpose() * l;
                ws = fabs(waveSpeed(ui, n)) * l;
                break;

            case Freestream:
            default:
                //Freestream Test
                block = flux(ui, c.uinf, n, Roe).transpose() * l;
                ws = fabs(max(waveSpeed(ui, n), waveSpeed(c.uinf, n))) * l;
                break;


            }
            //cout << block << endl;
#pragma omp critical 
            {
                R[j * u.u.cols() + 0] += block(0, 0);
                R[j * u.u.cols() + 1] += block(0, 1);
                R[j * u.u.cols() + 2] += block(0, 2);
                R[j * u.u.cols() + 3] += block(0, 3);

                sl[j] += ws;
            }
        }
    }

}


Matrix residual(FVstate& u, FVmesh& m, FVConditions& c, Matrix& dt) {
    //initialize

    double* Rval = new double[u.u.size()];
    double* sl = new double[u.u.rows()];
    int i, j;

    memset(Rval, 0, sizeof(double) * u.u.size());
    memset(sl, 0, sizeof(double) * u.u.rows());

    //gradient for second order
    if (u.Order == 2)
        gradient(u, m);

    //internal edges
    interioredges(u, m, Rval, sl);

    //boundary edges
    boundaryedges(u, m, c, Rval, sl);

    Matrix R(Rval, u.u.rows(), u.u.cols());

    //set time step
    dt = Matrix(u.u.rows(), 4);
    for (i = 0; i < u.u.rows(); ++i) {
        dt(i, 0) = 2 * CFL / sl[i];
        for (j = 1; j < dt.cols(); ++j) {
            dt(i, j) = dt(i, 0);

        }
    }
    //dt.print();
    delete[] Rval;
    delete[] sl;
    return R;
}


void grad_internal(FVstate& u, FVmesh& m, double* gradux, double* graduy) {

    Matrix goin_x = Matrix(1, 4);
    Matrix goin_y = Matrix(1, 4);
    Matrix stateAdd = Matrix(1, 4);
    Matrix ui = Matrix(1, 4);
    Matrix n = Matrix(1, 2);
    int i, j, k;

#pragma omp parallel 
    {
#pragma omp for private(j,k,stateAdd,n,goin_x,goin_y)
        for (i = 0; i < m.I2E.rows(); ++i) {
            j = m.I2E(i, 0);
            k = m.I2E(i, 2);

            stateAdd = u.u.getBlock(j, 0, 1, 4) + u.u.getBlock(k, 0, 1, 4);

            n = m.In.getBlock(i, 0, 1, 2) * (m.Il(i, 0) / 2);

            goin_x = stateAdd * n(0, 0);
            goin_y = stateAdd * n(0, 1);

#pragma omp critical
            {
                gradux[j * u.u.cols() + 0] += goin_x(0, 0);
                gradux[j * u.u.cols() + 1] += goin_x(0, 1);
                gradux[j * u.u.cols() + 2] += goin_x(0, 2);
                gradux[j * u.u.cols() + 3] += goin_x(0, 3);

                gradux[k * u.u.cols() + 0] -= goin_x(0, 0);
                gradux[k * u.u.cols() + 1] -= goin_x(0, 1);
                gradux[k * u.u.cols() + 2] -= goin_x(0, 2);
                gradux[k * u.u.cols() + 3] -= goin_x(0, 3);


                graduy[j * u.u.cols() + 0] += goin_y(0, 0);
                graduy[j * u.u.cols() + 1] += goin_y(0, 1);
                graduy[j * u.u.cols() + 2] += goin_y(0, 2);
                graduy[j * u.u.cols() + 3] += goin_y(0, 3);

                graduy[k * u.u.cols() + 0] -= goin_y(0, 0);
                graduy[k * u.u.cols() + 1] -= goin_y(0, 1);
                graduy[k * u.u.cols() + 2] -= goin_y(0, 2);
                graduy[k * u.u.cols() + 3] -= goin_y(0, 3);
            }
        }
    }
}

void grad_boundary(FVstate& u, FVmesh& m, double* gradux, double* graduy) {


    Matrix goin_x = Matrix(1, 4);
    Matrix goin_y = Matrix(1, 4);
    Matrix ui = Matrix(1, 4);
    Matrix n = Matrix(1, 2);
    int i, j;

#pragma omp parallel 
    {
#pragma omp for private(j,ui,n,goin_x,goin_y)
        for (i = 0; i < m.B2E.rows(); ++i) {
            j = m.B2E(i, 0);
            ui = u.u.getBlock(j, 0, 1, 4);
            n = m.Bn.getBlock(i, 0, 1, 2) * m.Bl(i, 0);

            goin_x = ui * n(0, 0);
            goin_y = ui * n(0, 1);

#pragma omp critical
            {
                gradux[j * u.u.cols() + 0] += goin_x(0, 0);
                gradux[j * u.u.cols() + 1] += goin_x(0, 1);
                gradux[j * u.u.cols() + 2] += goin_x(0, 2);
                gradux[j * u.u.cols() + 3] += goin_x(0, 3);


                graduy[j * u.u.cols() + 0] += goin_y(0, 0);
                graduy[j * u.u.cols() + 1] += goin_y(0, 1);
                graduy[j * u.u.cols() + 2] += goin_y(0, 2);
                graduy[j * u.u.cols() + 3] += goin_y(0, 3);
            }
        }
    }

}

void gradient(FVstate& u, FVmesh& m) {
    //reset gradient of u

    double* gradux = new double[u.u.size()];
    double* graduy = new double[u.u.size()];

    memset(gradux, 0, sizeof(double) * u.u.size());
    memset(graduy, 0, sizeof(double) * u.u.size());

    //internal edges
    grad_internal(u, m, gradux, graduy);

    //boundary edges
    grad_boundary(u, m, gradux, graduy);

    u.gradux = Matrix(gradux, u.u.rows(), u.u.cols());
    u.graduy = Matrix(graduy, u.u.rows(), u.u.cols());

    u.gradux = u.gradux % m.A;
    u.graduy = u.graduy % m.A;

    delete[] gradux;
    delete[] graduy;

    return;
}
