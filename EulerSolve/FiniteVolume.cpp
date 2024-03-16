#include "Tools.h"
#include <omp.h>

Matrix FV_solve(FVstate& u, FVmesh m, FVConditions c) {
    //constants
    double CFL = 0.5;
    double tol = 1e-7;
    int MaxIter = 100000;

    //initialize
    Matrix e(MaxIter + 1, 1);
    Matrix dt(u.u.rows(), u.u.cols());
    Matrix R(u.u.rows(), u.u.cols());
    int iter = 0;

    //iterate
    if (u.Order == 1) {

        for (iter = 0; iter < MaxIter; ++iter) {
            R = residual(u, m, c, dt, CFL);
            u.u = u.u - dt % R;
            e(iter, 0) = R.max();
            cout << e(iter, 0) << endl;
            if (e(iter, 0) < tol)
                break;
            if (e(iter, 0) > 1e10 || isnan(e(iter, 0)))
                break;
            if (iter % 1000 == 0)
                printResults(u.u, e);
        }
    } else if (u.Order == 2) {

        FVstate ufe(2, m, c);
        Matrix dtdummy(u.u.rows(), u.u.cols());

        for (iter = 0; iter < MaxIter; ++iter) {
            R = residual(u, m, c, dt, CFL);
            ufe.u = u.u - dt % R;
            R = residual(ufe, m, c, dtdummy, CFL);
            u.u = (u.u + ufe.u - dt % R) * .5;
            e(iter, 0) = R.max();
            cout << e(iter, 0) << endl;
            if (e(iter, 0) < tol)
                break;
            if (e(iter, 0) > 1e10 || isnan(e(iter, 0)))
                break;
        }
    }
    cout << "Stopped at iteration " << iter << endl;
    Matrix retVal(int(iter) + 1, 1);
    retVal = e.getBlock(0, 0, iter + 1, 1);
    return retVal;
}

Matrix residual(FVstate u, FVmesh m, FVConditions c, Matrix& dt, double CFL) {
    //initialize
    Matrix R(u.u.rows(), u.u.cols());
    Matrix sl(u.u.rows(), 1);


    Matrix gradux_l = Matrix(1, 4);
    Matrix graduy_l = Matrix(1, 4);
    Matrix gradux_r = Matrix(1, 4);
    Matrix graduy_r = Matrix(1, 4);
    Matrix gradux_i = Matrix(1, 4);
    Matrix graduy_i = Matrix(1, 4);
    Matrix uL       = Matrix(1, 4);
    Matrix uR       = Matrix(1, 4);
    Matrix ui       = Matrix(1, 4);
    Matrix uout     = Matrix(1, 4);
    Matrix block    = Matrix(1, 4);
    Matrix R_l      = Matrix(1, 4);
    Matrix R_r      = Matrix(1, 4);
    Matrix f        = Matrix(4, 1);
    Matrix n        = Matrix(1, 2);
    Matrix uin      = Matrix(1, 2);
    int i, j, k;
    double l;

    //gradient for second order
    if (u.Order == 2)
        gradient(u, m);

    //internal edges
#pragma omp parallel
    {
#pragma omp for reduction(+:R,sl)
        for (i = 0; i < m.I2E.rows(); ++i) {


            j = int(m.I2E(i, 0));
            k = int(m.I2E(i, 2));
            uL = u.u.getBlock(j, 0, 1, 4);
            uR = u.u.getBlock(k, 0, 1, 4);

            if (u.Order == 2) {
                gradux_l = u.gradux.getBlock(j, 0, 1, 4);
                graduy_l = u.graduy.getBlock(j, 0, 1, 4);
                gradux_r = u.gradux.getBlock(k, 0, 1, 4);
                graduy_r = u.graduy.getBlock(k, 0, 1, 4);

                //std::cout << "before\n" << uL << endl << uR << endl;

                uL += (gradux_l * (m.Ir(i, 0) - m.C(j, 0)) + graduy_l * (m.Ir(i, 1) - m.C(j, 1)));
                uR += (gradux_r * (m.Ir(i, 0) - m.C(k, 0)) + graduy_r * (m.Ir(i, 1) - m.C(k, 1)));


                //std::cout << "after\n" << uL << endl << uR << endl;

                //std::cout << gradux_l << endl;
                //std::cout << graduy_l << endl;
                //std::cout << gradux_r << endl;
                //std::cout << graduy_r << endl;
                //std::cout << endl;


            }
            n = m.In.getBlock(i, 0, 1, 2);
            l = m.Il(i, 0);
            f = flux(uL, uR, n, 0);
            //f.print();
            double ws = max(waveSpeed(uL, n), waveSpeed(uR, n));

            R_l = R.getBlock(j, 0, 1, 4) + (f.transpose() * l);
            R_r = R.getBlock(k, 0, 1, 4) - (f.transpose() * l);

            R.setBlock(j, 0, R_l);
            R.setBlock(k, 0, R_r);

            //vec d = R.getBlock(j, 0, 1, 4);
            //d.print();

            sl(j, 0) = sl(j, 0) + ws * l;
            sl(k, 0) = sl(k, 0) + ws * l;

        }
    }
    //R.print();

    //boundary edges
    for (i = 0; i < m.B2E.rows(); ++i) {
        j = int(m.B2E(i, 0));
        ui = u.u.getBlock(j, 0, 1, 4);
        if (u.Order == 2) {
            gradux_i = u.gradux.getBlock(j, 0, 1, 4);
            graduy_i = u.graduy.getBlock(j, 0, 1, 4);
            ui += (gradux_i * (m.Br(i, 0) - m.C(j, 0)) + graduy_i * (m.Br(i, 1) - m.C(j, 1)));
        }
        n = m.Bn.getBlock(i, 0, 1, 2);

        l = m.Bl(i, 0);
        uout(1, 4);
        uin(1, 4);

        block = R.getBlock(j, 0, 1, 4);
        switch ((int)m.B2E(i, 2)) {
        case 1:
        case 3:
            //Inviscid Wall
            //cout << j << endl;
            //block.print();
            block += wallFlux(ui, n) * l;
            R.setBlock(j, 0, block);
            sl(j, 0) = sl(j, 0) + fabs(waveSpeed(ui, n)) * l;
            //cout << endl << "wall" << endl;
            //R.print();
            break;

        case 2:
            //Subsonic Outlet
            uout = Outflow(ui, n, c.Pinf);
            block += (F(uout) * (n.transpose())).transpose() * l;
            R.setBlock(j, 0, block);
            sl(j, 0) = sl(j, 0) + fabs(waveSpeed(uout, n)) * l;
            //cout << endl << "out" << endl;
            //R.print();
            break;

        case 4:
            //Subsonic Inlet
            //cout << j << endl;
            uin = Inflow(ui, n, c.Tt, c.Pt, c.a, c.R);
            block += (F(uin) * (n.transpose())).transpose() * l;
            R.setBlock(j, 0, block);
            sl(j, 0) = sl(j, 0) + fabs(waveSpeed(uin, n)) * l;
            //cout << endl << "in" << endl;
            //R.print();
            break;

        case 5:
            //Supersonic Inlet
            //block.print();
            //cout << j << endl;
            block += flux(ui, c.uinf, n, 0).transpose() * l;
            R.setBlock(j, 0, block);
            sl(j, 0) = sl(j, 0) + fabs(max(waveSpeed(ui, n), waveSpeed(c.uinf, n))) * l;
            break;

        case 6:
            //Supersonic Outlet
            //block.print();
            //cout << j << endl;
            block += (F(ui) * (n.transpose())).transpose() * l;
            R.setBlock(j, 0, block);
            sl(j, 0) = sl(j, 0) + fabs(waveSpeed(ui, n)) * l;
            break;

        default:
            //Freestream Test
            block += flux(ui, c.uinf, n, 0).transpose() * l;
            R.setBlock(j, 0, block);
            sl(j, 0) = sl(j, 0) + fabs(max(waveSpeed(ui, n), waveSpeed(c.uinf, n))) * l;
            break;

        }
    }
    //R.print();
    //set time step
    dt = Matrix(sl.rows(), 4);
    for (i = 0; i < sl.rows(); ++i) {
        dt(i, 0) = 2 * CFL / sl(i, 0);
        for (j = 1; j < dt.cols(); ++j) {
            dt(i, j) = dt(i, 0);

        }
    }
    //dt.print();
    return R;
}

void gradient(FVstate& u, FVmesh& m) {
    //reset gradient of u
    u.gradux = Matrix(u.u.rows(), u.u.cols());
    u.graduy = Matrix(u.u.rows(), u.u.cols());

    Matrix goin_lx = Matrix(1, 4);
    Matrix goin_rx = Matrix(1, 4);
    Matrix goin_ly = Matrix(1, 4);
    Matrix goin_ry = Matrix(1, 4);
    Matrix goin_x = Matrix(1, 4);
    Matrix goin_y = Matrix(1, 4);
    Matrix stateAdd = Matrix(1, 4);
    Matrix ui = Matrix(1, 4);
    Matrix n = Matrix(1, 2);

    //internal edges
    for (int i = 0; i < m.I2E.rows(); ++i) {
        int j = int(m.I2E(i, 0));
        int k = int(m.I2E(i, 2));
        //Matrix uL = u.u.getBlock(j, 0, 1, 4);
        //Matrix uR = u.u.getBlock(k, 0, 1, 4);
        stateAdd = u.u.getBlock(j, 0, 1, 4) + u.u.getBlock(k, 0, 1, 4);
        
        n = m.In.getBlock(i, 0, 1, 2) * (m.Il(i, 0) / 2);
        //double l = m.Il(i, 0);


        //Matrix gradux_l = u.gradux.getBlock(j, 0, 1, 4);
        //Matrix graduy_l = u.graduy.getBlock(j, 0, 1, 4);
        //Matrix gradux_r = u.gradux.getBlock(k, 0, 1, 4);
        //Matrix graduy_r = u.graduy.getBlock(k, 0, 1, 4);

        goin_lx = u.gradux.getBlock(j, 0, 1, 4) + stateAdd * n(0, 0);
        goin_rx = u.gradux.getBlock(k, 0, 1, 4) - stateAdd * n(0, 0);
        goin_ly = u.graduy.getBlock(j, 0, 1, 4) + stateAdd * n(0, 1);
        goin_ry = u.graduy.getBlock(k, 0, 1, 4) - stateAdd * n(0, 1);

        u.gradux.setBlock(j, 0, goin_lx);
        u.gradux.setBlock(k, 0, goin_rx);

        u.graduy.setBlock(j, 0, goin_ly);
        u.graduy.setBlock(k, 0, goin_ry);
    }

    //boundary edges
    for (int i = 0; i < m.B2E.rows(); ++i) {
        int j = int(m.B2E(i, 0));
        ui = u.u.getBlock(j, 0, 1, 4);
        n = m.Bn.getBlock(i, 0, 1, 2) * m.Bl(i, 0);
        //double l = m.Bl(i, 0);

        //Matrix gradux_i = u.gradux.getBlock(j, 0, 1, 4);
        //Matrix graduy_i = u.graduy.getBlock(j, 0, 1, 4);

        goin_x = u.gradux.getBlock(j, 0, 1, 4) + ui * n(0, 0);
        goin_y = u.graduy.getBlock(j, 0, 1, 4) + ui * n(0, 1);

        u.gradux.setBlock(j, 0, goin_x);
        u.graduy.setBlock(j, 0, goin_y);

    }

    //std::cout << u.gradux << endl;
    //std::cout << u.graduy << endl;
    //std::cout << A << endl;

    u.gradux = u.gradux % m.A;
    u.graduy = u.graduy % m.A;

    return;
}
