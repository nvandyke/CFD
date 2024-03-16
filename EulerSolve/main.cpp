#include "Tools.h"
#include "matrix.h"
#include <iostream>
#include <omp.h>
#include <chrono>


void FVrun() {
    //string filename = "mesh.txt";
    FVmesh m("meshes/bump0.gri");
    double Mach = 0.5;
    double angleOfAttack = 0.0;
    int order = 1;

    FVConditions c(Mach, angleOfAttack);
    FVstate u(order, m, c);

    Matrix e = FV_solve(u, m, c);
    printResults(u.u, e);


    //u.Order = 2;
    //Matrix e2 = FV_solve(u, m, c);
    //printResults(u.u, e2);
    return;
}

//void FErun()
//{
//	string filename = "mesh.txt";
//	FEmesh m(filename);
//	FEConditions c(0.5, 0);
//	FEstate u(1, 1, m, c);
//	cout << m.nodes << endl;
//	return;
//}


int main() {
    omp_set_num_threads(8);
    std::cout.precision(15);
    std::cout << std::fixed;

    auto start_time = chrono::high_resolution_clock::now();
    FVrun();
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    std::cout << "took " << duration.count() << "s" << std::endl;

    std::cin.get();
    return 0;
}