#include "Tools.h"
#include <omp.h>


extern "C" {
    void FVrun(std::string, std::string, int, double, double);
    void cudaStart(int);
    void cudaEnd();
}


void FVrun(std::string ic, std::string mesh, int order, double Mach, double angleOfAttack) {
    //cudaStart(9);

    //cudaEnd();
    //return;
    Mesh m(mesh);

    State u(order, m, Mach, angleOfAttack);

    if (ic != "") {
        std::ifstream inputfile;
        inputfile.open(ic);
        inputfile >> u.u;
        inputfile.close();
    }

    Matrix e = u.FV_solve();
    printResults(u.u, e);
    return;
}


Matrix State::FV_solve() {
    //initialize
    cudaStart(u.size());
    Matrix e(MaxIter + 1, 1);
    Matrix dt(u.rows(), u.cols());
    Matrix R(u.rows(), u.cols());
    int iter = 0;

    //higher-order initialize
    State ufe(Order, mesh, M, a);
    Matrix dtdummy(u.rows(), u.cols());
    Matrix R2(u.rows(), u.cols());
    Matrix R3(u.rows(), u.cols());
    Matrix R4(u.rows(), u.cols());
    Matrix bigR(u.rows(), u.cols());

    //iterate
    for (iter = 0; iter < MaxIter; ++iter) {
        R = residual(dt);
        
        if (Order == 1) {
            u = u - dt % R;
            //u = u - dt.multInPlace(R);
        } else if (Order == 2) {
            ufe.u = u - dt % R;
            R2 = ufe.residual(dtdummy);
            u = (u + ufe.u - dt % R2) * .5;
        } else if (Order < 5) {
            ufe.u = u - .5 * (dt % R);
            R2 = ufe.residual(dtdummy);
            ufe.u = u - .5 * (dt % R2);
            R3 = ufe.residual(dtdummy);
            ufe.u = u - dt % R3;
            R4 = ufe.residual(dtdummy);
            bigR = R + 2 * (R2 + R3) + R4;
            u = u - (dt % bigR) / 6;
        }

        e(iter, 0) = std::max(R.max(), -R.min());
        std::cout << e(iter, 0) << std::endl;
        if (e(iter, 0) < tol)
            break;
        if (e(iter, 0) > max || isnan(e(iter, 0)))
            break;
        if (iter % 1000 == 0)
            printResults(u, e);
    }

    cudaEnd();
    std::cout << "Stopped at iteration " << iter << std::endl;
    Matrix retVal(int(iter) + 1, 1);
    retVal = e.getBlock(0, 0, iter + 1, 1);
    return retVal;
}

void State::interioredges(double* R, double* sl) {


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
        for (i = 0; i < mesh.I2E.rows(); ++i) {

            j = mesh.I2E(i, 0);
            k = mesh.I2E(i, 2);
            uL = u.getBlock(j, 0, 1, 4);
            uR = u.getBlock(k, 0, 1, 4);

            if (Order == 2) {
                gradux_l = gradux.getBlock(j, 0, 1, 4);
                graduy_l = graduy.getBlock(j, 0, 1, 4);
                gradux_r = gradux.getBlock(k, 0, 1, 4);
                graduy_r = graduy.getBlock(k, 0, 1, 4);

                uL += (gradux_l * (mesh.Ir(i, 0) - mesh.C(j, 0)) + graduy_l * (mesh.Ir(i, 1) - mesh.C(j, 1)));
                uR += (gradux_r * (mesh.Ir(i, 0) - mesh.C(k, 0)) + graduy_r * (mesh.Ir(i, 1) - mesh.C(k, 1)));

            }
            n = mesh.In.getBlock(i, 0, 1, 2);
            l = mesh.Il(i, 0);
            f = flux(uL, uR, n, flux_val);

            ws = std::max(waveSpeed(uL, n), waveSpeed(uR, n)) * l;

            fTl = f.transpose() * l;
#pragma omp critical 
            {
                R[u.cols() * j + 0] += fTl(0, 0);
                R[u.cols() * j + 1] += fTl(0, 1);
                R[u.cols() * j + 2] += fTl(0, 2);
                R[u.cols() * j + 3] += fTl(0, 3);

                R[u.cols() * k + 0] -= fTl(0, 0);
                R[u.cols() * k + 1] -= fTl(0, 1);
                R[u.cols() * k + 2] -= fTl(0, 2);
                R[u.cols() * k + 3] -= fTl(0, 3);

                sl[j] += ws;
                sl[k] += ws;
            }
        }
    }
}

void State::boundaryedges(double* R, double* sl) {


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
        for (i = 0; i < mesh.B2E.rows(); ++i) {
            j = mesh.B2E(i, 0);
            ui = u.getBlock(j, 0, 1, 4);
            if (Order == 2) {
                gradux_i = gradux.getBlock(j, 0, 1, 4);
                graduy_i = graduy.getBlock(j, 0, 1, 4);
                ui += (gradux_i * (mesh.Br(i, 0) - mesh.C(j, 0)) + graduy_i * (mesh.Br(i, 1) - mesh.C(j, 1)));
            }
            n = mesh.Bn.getBlock(i, 0, 1, 2);
            l = mesh.Bl(i, 0);

            switch (mesh.B2E(i, 2)) {
            case Inviscid_Wall:
                //Inviscid Wall
                block = wallFlux(ui, n) * l;
                ws = fabs(waveSpeed(ui, n)) * l;
                break;

            case Subsonic_Outlet:
                //Subsonic Outlet
                uout = Outflow(ui, n, Pinf);
                block = (F(uout) * (n.transpose())).transpose() * l;
                ws = fabs(waveSpeed(uout, n)) * l;
                break;

            case Subsonic_Inlet:
                //Subsonic Inlet
                uin = Inflow(ui, n, Tt, Pt, a, Rinf);
                block = (F(uin) * (n.transpose())).transpose() * l;
                ws = fabs(waveSpeed(uin, n)) * l;
                break;

            case Supersonic_Inlet:
                //Supersonic Inlet
                block = flux(ui, uinf, n, flux_val).transpose() * l;
                ws = fabs(std::max(waveSpeed(ui, n), waveSpeed(uinf, n))) * l;
                break;

            case Supersonic_Outlet:
                //Supersonic Outlet
                block = (F(ui) * (n.transpose())).transpose() * l;
                ws = fabs(waveSpeed(ui, n)) * l;
                break;

            case Freestream:
            default:
                //Freestream Test
                block = flux(ui, uinf, n, flux_val).transpose() * l;
                ws = fabs(std::max(waveSpeed(ui, n), waveSpeed(uinf, n))) * l;
                break;


            }
            //cout << block << endl;
#pragma omp critical 
            {
                R[j * u.cols() + 0] += block(0, 0);
                R[j * u.cols() + 1] += block(0, 1);
                R[j * u.cols() + 2] += block(0, 2);
                R[j * u.cols() + 3] += block(0, 3);

                sl[j] += ws;
            }
        }
    }

}


Matrix State::residual(Matrix& dt) {
    //initialize

    double* Rval = new double[u.size()];
    double* sl = new double[u.rows()];
    int i, j;

    memset(Rval, 0, sizeof(double) * u.size());
    memset(sl, 0, sizeof(double) * u.rows());

    //gradient for second order
    if (Order == 2)
        gradient();

    //internal edges
    interioredges(Rval, sl);

    //boundary edges
    boundaryedges(Rval, sl);

    Matrix R(Rval, u.rows(), u.cols());

    //set time step
    dt = Matrix(u.rows(), 4);
    for (i = 0; i < u.rows(); ++i) {
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


void State::grad_internal(double* gradux_l, double* graduy_l) {

    Matrix goin_x = Matrix(1, 4);
    Matrix goin_y = Matrix(1, 4);
    Matrix stateAdd = Matrix(1, 4);
    Matrix ui = Matrix(1, 4);
    Matrix n = Matrix(1, 2);
    int i, j, k;

#pragma omp parallel 
    {
#pragma omp for private(j,k,stateAdd,n,goin_x,goin_y)
        for (i = 0; i < mesh.I2E.rows(); ++i) {
            j = mesh.I2E(i, 0);
            k = mesh.I2E(i, 2);

            stateAdd = u.getBlock(j, 0, 1, 4) + u.getBlock(k, 0, 1, 4);

            n = mesh.In.getBlock(i, 0, 1, 2) * (mesh.Il(i, 0) / 2);

            goin_x = stateAdd * n(0, 0);
            goin_y = stateAdd * n(0, 1);

#pragma omp critical
            {
                gradux_l[j * u.cols() + 0] += goin_x(0, 0);
                gradux_l[j * u.cols() + 1] += goin_x(0, 1);
                gradux_l[j * u.cols() + 2] += goin_x(0, 2);
                gradux_l[j * u.cols() + 3] += goin_x(0, 3);

                gradux_l[k * u.cols() + 0] -= goin_x(0, 0);
                gradux_l[k * u.cols() + 1] -= goin_x(0, 1);
                gradux_l[k * u.cols() + 2] -= goin_x(0, 2);
                gradux_l[k * u.cols() + 3] -= goin_x(0, 3);


                graduy_l[j * u.cols() + 0] += goin_y(0, 0);
                graduy_l[j * u.cols() + 1] += goin_y(0, 1);
                graduy_l[j * u.cols() + 2] += goin_y(0, 2);
                graduy_l[j * u.cols() + 3] += goin_y(0, 3);

                graduy_l[k * u.cols() + 0] -= goin_y(0, 0);
                graduy_l[k * u.cols() + 1] -= goin_y(0, 1);
                graduy_l[k * u.cols() + 2] -= goin_y(0, 2);
                graduy_l[k * u.cols() + 3] -= goin_y(0, 3);
            }
        }
    }
}

void State::grad_boundary(double* gradux_l, double* graduy_l) {


    Matrix goin_x = Matrix(1, 4);
    Matrix goin_y = Matrix(1, 4);
    Matrix ui = Matrix(1, 4);
    Matrix n = Matrix(1, 2);
    int i, j;

#pragma omp parallel 
    {
#pragma omp for private(j,ui,n,goin_x,goin_y)
        for (i = 0; i < mesh.B2E.rows(); ++i) {
            j = mesh.B2E(i, 0);
            ui = u.getBlock(j, 0, 1, 4);
            n = mesh.Bn.getBlock(i, 0, 1, 2) * mesh.Bl(i, 0);

            goin_x = ui * n(0, 0);
            goin_y = ui * n(0, 1);

#pragma omp critical
            {
                gradux_l[j * u.cols() + 0] += goin_x(0, 0);
                gradux_l[j * u.cols() + 1] += goin_x(0, 1);
                gradux_l[j * u.cols() + 2] += goin_x(0, 2);
                gradux_l[j * u.cols() + 3] += goin_x(0, 3);


                graduy_l[j * u.cols() + 0] += goin_y(0, 0);
                graduy_l[j * u.cols() + 1] += goin_y(0, 1);
                graduy_l[j * u.cols() + 2] += goin_y(0, 2);
                graduy_l[j * u.cols() + 3] += goin_y(0, 3);
            }
        }
    }

}

void State::gradient() {
    //reset gradient of u

    double* gradux_local = new double[u.size()];
    double* graduy_local = new double[u.size()];

    memset(gradux_local, 0, sizeof(double) * u.size());
    memset(graduy_local, 0, sizeof(double) * u.size());

    //internal edges
    grad_internal(gradux_local, graduy_local);

    //boundary edges
    grad_boundary(gradux_local, graduy_local);

    gradux = Matrix(gradux_local, u.rows(), u.cols());
    graduy = Matrix(graduy_local, u.rows(), u.cols());

    gradux = gradux % mesh.A;
    graduy = graduy % mesh.A;

    delete[] gradux_local;
    delete[] graduy_local;

    return;
}
