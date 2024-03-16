#include "Tools.h"
#include <omp.h>
/*
VectorXd FV_solve(FVstate& u, const FVmesh& m, const FVConditions& c)
{
    //constants
    double CFL = 0.5;
    double tol = 1e-7;
    int MaxIter = 500000;

    //initialize
    VectorXd e(MaxIter);
    MatrixXd dt(u.u.rows(), u.u.cols());
    MatrixXd R(u.u.rows(), u.u.cols());
    int iter = 0;

    //iterate
    if (u.Order == 1) {

        for (iter = 0; iter < MaxIter; ++iter) {
            R = residual(u, m, c, dt, CFL);
            u.u.noalias() -= dt.cwiseProduct(R);
            e(iter) = max(R.maxCoeff(), (-R).maxCoeff());
            cout << e(iter) << endl;
            if (e(iter) < tol)
                break;
            if (e(iter) > 1e10 || isnan(e(iter)))
                break;
        }
    }
    else if (u.Order == 2) {

        FVstate ufe(2, m, c);
        MatrixXd dtdummy(u.u.rows(), u.u.cols());

        for (iter = 0; iter < MaxIter; ++iter) {
            R = residual(u, m, c, dt, CFL);
            ufe.u.noalias() = u.u - dt.cwiseProduct(R);
            R = residual(ufe, m, c, dtdummy, CFL);
            u.u.noalias() = .5 * (u.u + ufe.u - dt.cwiseProduct(R)).eval();
            e(iter) = max(R.maxCoeff(), (-R).maxCoeff());
            cout << e(iter) << endl;
            if (e(iter) < tol)
                break;
            if (e(iter) > 1e10 || isnan(e(iter)))
                break;
        }
    }
    cout << "Stopped at iteration " << iter << endl;
    return e.head(iter + 1);
}
*/

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
/*
MatrixXd residual(FVstate& u, const FVmesh& m, const FVConditions& c, MatrixXd& dt, double CFL)
{
    //initialize
    MatrixXd R;
    R.setZero(u.u.rows(), u.u.cols());
    dt = R;
    MatrixXd sl;
    sl.setZero(u.u.rows(), 1);

    //gradient for second order
    if (u.Order == 2)
        gradient(u, m);

    //internal edges
    for (int i = 0; i < m.I2E.rows(); ++i) {
        int j = m.I2E(i, 0) - 1;
        int k = m.I2E(i, 2) - 1;
        Vector4d uL = u.u.block<1, 4>(j, 0);
        Vector4d uR = u.u.block<1, 4>(k, 0);
        if (u.Order == 2) {
            uL += u.gradux.block<1, 4>(j, 0) * (m.Ir(i, 0) - m.C(j, 0)) + u.graduy.block<1, 4>(j, 0) * (m.Ir(i, 1) - m.C(j, 1));
            uR += u.gradux.block<1, 4>(k, 0) * (m.Ir(i, 0) - m.C(k, 0)) + u.graduy.block<1, 4>(k, 0) * (m.Ir(i, 1) - m.C(k, 1));
        }
        Vector2d n(m.In(i, 0), m.In(i, 1));
        double l = m.Il(i);
        Vector4d f = flux(uL, uR, n, 0);
        double ws = fabs(max(waveSpeed(uL, n), waveSpeed(uR, n)));

        R.block<1, 4>(j, 0).noalias() += f * l;
        R.block<1, 4>(k, 0).noalias() -= f * l;

        sl(j) += ws * l;
        sl(k) += ws * l;
    }

    //boundary edges
    for (int i = 0; i < m.B2E.rows(); ++i) {
        int j = m.B2E(i, 0) - 1;
        Vector4d ui = u.u.block<1, 4>(j, 0);
        if (u.Order == 2) {
            ui += u.gradux.block<1, 4>(j, 0) * (m.Br(i, 0) - m.C(j, 0)) + u.graduy.block<1, 4>(j, 0) * (m.Br(i, 1) - m.C(j, 1));
        }
        Vector2d n(m.Bn(i, 0), m.Bn(i, 1));
        double l = m.Bl(i);
        Vector4d uout(0, 0, 0, 0);
        Vector4d uin;

        switch (m.B2E(i, 2)) {
        case 1:
        case 3:
            //Inviscid Wall
            R.block<1, 4>(j, 0).noalias() += wallFlux(ui, n) * l;
            sl(j) += fabs(waveSpeed(ui, n)) * l;
            break;

        case 2:
            //Subsonic Outlet
            uout = Outflow(ui, n, c.Pinf);
            R.block<1, 4>(j, 0).noalias() += F(uout) * n * l;
            sl(j) += fabs(waveSpeed(uout, n)) * l;
            break;

        case 4:
            //Subsonic Inlet
            uin = Inflow(ui, n, c.Tt, c.Pt, c.a, c.R);
            R.block<1, 4>(j, 0).noalias() += F(uin) * n * l;
            sl(j) += fabs(waveSpeed(uin, n)) * l;
            break;

        case 5:
            //Supersonic Inlet
            R.block<1, 4>(j, 0).noalias() += flux(ui, c.uinf, n, 0) * l;
            sl(j) += fabs(max(waveSpeed(ui, n), waveSpeed(c.uinf, n))) * l;
            break;

        case 6:
            //Supersonic Outlet
            R.block<1, 4>(j, 0).noalias() += F(ui) * n * l;
            sl(j) += fabs(waveSpeed(ui, n)) * l;
            break;

        default:
            //Freestream Test
            R.block<1, 4>(j, 0).noalias() += flux(ui, c.uinf, n, 0) * l;
            sl(j) += fabs(max(waveSpeed(ui, n), waveSpeed(c.uinf, n))) * l;
            break;

        }
    }

    //set time step
    for (int i = 0; i < sl.rows(); ++i) {
        dt(i, 0) = 2 * CFL / sl(i);
        for (int j = 1; j < 4; ++j) {
            dt(i, j) = dt(i, 0);
        }
    }

    return R;
}
*/
Matrix residual(FVstate u, FVmesh m, FVConditions c, Matrix& dt, double CFL) {
    //initialize
    Matrix R(u.u.rows(), u.u.cols());
    Matrix sl(u.u.rows(), 1);

    //gradient for second order
    if (u.Order == 2)
        gradient(u, m);

    //internal edges
    for (int i = 0; i < m.I2E.rows(); ++i) {


        int j = int(m.I2E(i, 0));
        int k = int(m.I2E(i, 2));
        Matrix uL = u.u.getBlock(j, 0, 1, 4);
        Matrix uR = u.u.getBlock(k, 0, 1, 4);

        if (u.Order == 2) {
            Matrix gradux_l = u.gradux.getBlock(j, 0, 1, 4);
            Matrix graduy_l = u.graduy.getBlock(j, 0, 1, 4);
            Matrix gradux_r = u.gradux.getBlock(k, 0, 1, 4);
            Matrix graduy_r = u.graduy.getBlock(k, 0, 1, 4);

            //std::cout << "before\n" << uL << endl << uR << endl;

            uL = uL + (gradux_l * (m.Ir(i, 0) - m.C(j, 0)) + graduy_l * (m.Ir(i, 1) - m.C(j, 1)));
            uR = uR + (gradux_r * (m.Ir(i, 0) - m.C(k, 0)) + graduy_r * (m.Ir(i, 1) - m.C(k, 1)));


            //std::cout << "after\n" << uL << endl << uR << endl;

            //std::cout << gradux_l << endl;
            //std::cout << graduy_l << endl;
            //std::cout << gradux_r << endl;
            //std::cout << graduy_r << endl;
            //std::cout << endl;


        }
        Matrix n = m.In.getBlock(i, 0, 1, 2);
        double l = m.Il(i, 0);
        Matrix f = flux(uL, uR, n, 0);
        //f.print();
        double ws = max(waveSpeed(uL, n), waveSpeed(uR, n));

        Matrix R_l = R.getBlock(j, 0, 1, 4) + (f.transpose() * l);
        Matrix R_r = R.getBlock(k, 0, 1, 4) - (f.transpose() * l);

        R.setBlock(j, 0, R_l);
        R.setBlock(k, 0, R_r);

        //vec d = R.getBlock(j, 0, 1, 4);
        //d.print();

        sl(j, 0) = sl(j, 0) + ws * l;
        sl(k, 0) = sl(k, 0) + ws * l;

    }
    //R.print();

    //boundary edges
    for (int i = 0; i < m.B2E.rows(); ++i) {
        int j = int(m.B2E(i, 0));
        Matrix ui = u.u.getBlock(j, 0, 1, 4);
        if (u.Order == 2) {
            Matrix gradux_i = u.gradux.getBlock(j, 0, 1, 4);
            Matrix graduy_i = u.graduy.getBlock(j, 0, 1, 4);
            ui = ui + (gradux_i * (m.Br(i, 0) - m.C(j, 0)) + graduy_i * (m.Br(i, 1) - m.C(j, 1)));
        }
        Matrix n = m.Bn.getBlock(i, 0, 1, 2);
        //n = n.transpose();
        double l = m.Bl(i, 0);
        Matrix uout(1, 4);
        Matrix uin(1, 4);

        Matrix block = R.getBlock(j, 0, 1, 4);
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
    for (int i = 0; i < sl.rows(); ++i) {
        dt(i, 0) = 2 * CFL / sl(i, 0);
        //if (dt.getAt(i, 0) > 10) {
        //	cout << i << endl;
        //}
        for (int j = 1; j < 4; ++j) {
            dt(i, j) = dt(i, 0);

        }
    }
    //dt.print();
    return R;
}
/*
void gradient(FVstate& u, FVmesh& m)
{
    //reset gradient of u
    u.gradux.setZero(u.u.rows(), u.u.cols());
    u.graduy.setZero(u.u.rows(), u.u.cols());

    //internal edges
    for (int i = 0; i < m.I2E.rows(); ++i) {
        int j = m.I2E(i, 0) - 1;
        int k = m.I2E(i, 2) - 1;
        Vector4d uL = u.u.block<1, 4>(j, 0);
        Vector4d uR = u.u.block<1, 4>(k, 0);
        Vector2d n(m.In(i, 0), m.In(i, 1));
        double l = m.Il(i);

        u.gradux.block<1, 4>(j, 0).noalias() += .5 * (uL + uR) * n(0) * l;
        u.gradux.block<1, 4>(k, 0).noalias() -= .5 * (uL + uR) * n(0) * l;

        u.graduy.block<1, 4>(j, 0).noalias() += .5 * (uL + uR) * n(1) * l;
        u.graduy.block<1, 4>(k, 0).noalias() -= .5 * (uL + uR) * n(1) * l;
    }

    //boundary edges
    for (int i = 0; i < m.B2E.rows(); ++i) {
        int j = m.B2E(i, 0) - 1;
        Vector4d ui = u.u.block<1, 4>(j, 0);
        Vector2d n(m.Bn(i, 0), m.Bn(i, 1));
        double l = m.Bl(i);

        u.gradux.block<1, 4>(j, 0).noalias() += ui * n(0) * l;
        u.graduy.block<1, 4>(j, 0).noalias() += ui * n(1) * l;
    }

    //normalize to area
    MatrixXd A(m.A.rows(), 4);
    for (int i = 0; i < A.rows(); ++i) {
        A(i, 0) = 1 / m.A(i);
        for (int j = 1; j < 4; ++j) {
            A(i, j) = A(i, 0);
        }
    }

    u.gradux.noalias() = u.gradux.cwiseProduct(A);
    u.graduy.noalias() = u.graduy.cwiseProduct(A);

    return;
}
*/
void gradient(FVstate& u, FVmesh& m) {
    //reset gradient of u
    u.gradux = Matrix(u.u.rows(), u.u.cols());
    u.graduy = Matrix(u.u.rows(), u.u.cols());

    //internal edges
    for (int i = 0; i < m.I2E.rows(); ++i) {
        int j = int(m.I2E(i, 0));
        int k = int(m.I2E(i, 2));
        Matrix uL = u.u.getBlock(j, 0, 1, 4);
        Matrix uR = u.u.getBlock(k, 0, 1, 4);
        Matrix n = m.In.getBlock(i, 0, 1, 2);
        double l = m.Il(i, 0);


        Matrix gradux_l = u.gradux.getBlock(j, 0, 1, 4);
        Matrix graduy_l = u.graduy.getBlock(j, 0, 1, 4);
        Matrix gradux_r = u.gradux.getBlock(k, 0, 1, 4);
        Matrix graduy_r = u.graduy.getBlock(k, 0, 1, 4);

        Matrix goin_lx = gradux_l + (uL + uR) * (.5 * n(0, 0) * l);
        Matrix goin_rx = gradux_r - (uL + uR) * (.5 * n(0, 0) * l);
        Matrix goin_ly = graduy_l + (uL + uR) * (.5 * n(0, 1) * l);
        Matrix goin_ry = graduy_r - (uL + uR) * (.5 * n(0, 1) * l);
        u.gradux.setBlock(j, 0, goin_lx);
        u.gradux.setBlock(k, 0, goin_rx);

        u.graduy.setBlock(j, 0, goin_ly);
        u.graduy.setBlock(k, 0, goin_ry);
    }

    //boundary edges
    for (int i = 0; i < m.B2E.rows(); ++i) {
        int j = int(m.B2E(i, 0));
        Matrix ui = u.u.getBlock(j, 0, 1, 4);
        Matrix n = m.Bn.getBlock(i, 0, 1, 2);
        double l = m.Bl(i, 0);

        Matrix gradux_i = u.gradux.getBlock(j, 0, 1, 4);
        Matrix graduy_i = u.graduy.getBlock(j, 0, 1, 4);

        Matrix goin_x = gradux_i + ui * n(0, 0) * l;
        Matrix goin_y = graduy_i + ui * n(0, 1) * l;
        u.gradux.setBlock(j, 0, goin_x);
        u.graduy.setBlock(j, 0, goin_y);

    }

    //normalize to area
    Matrix A(m.A.rows(), 4);
    for (int i = 0; i < A.rows(); ++i) {
        A(i, 0) = 1 / m.A(i, 0);
        for (int j = 1; j < 4; ++j) {
            A(i, j) = A(i, 0);
        }
    }

    //std::cout << u.gradux << endl;
    //std::cout << u.graduy << endl;
    //std::cout << A << endl;

    u.gradux = u.gradux % A;
    u.graduy = u.graduy % A;

    return;
}
