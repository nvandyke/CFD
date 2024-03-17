#include "Tools.h"
#include "matrix.h"
#include <iostream>
#include <omp.h>
#include <chrono>
#include <string>


void FVrun(string ic, string mesh) {
    //string filename = "mesh.txt";
    FVmesh m(mesh);
    double Mach = 0.5;
    double angleOfAttack = 0.0;
    int order = 1;

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
inline bool exists(const string& name) {
    ifstream f(name.c_str());
    return f.good();
}


int main(int argc, char *argv[]) {
    omp_set_num_threads(8);
    std::cout.precision(15);
    std::cout << std::fixed;
    string initialfile = "";
    string meshfile = "meshes\\bump0.gri";

    while ((++argv)[0]) {
        if (argv[0][0] == '-') {
            switch (argv[0][1]) {

            
            case 'i':
                if (exists((++argv)[0])) {
                    initialfile = argv[0];
                    std::cout << "init file: " << initialfile << endl;
                }
                break;
            case 'm':
                if (exists((++argv)[0])) {
                    meshfile = argv[0];
                    std::cout << "mesh file: " << meshfile << endl;
                }
                break;
            case 'h':
            default:
                cout << "options are:" << endl;
                cout << "\t-i: input file" << endl;
                cout << "\t-m: mesh file" << endl;
                cout << "\t-h: help menu" << endl;
                return 1;
                break;
            }
        }

    }

    auto start_time = chrono::high_resolution_clock::now();
    FVrun(initialfile,meshfile);
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    std::cout << "took " << duration.count() << "s" << std::endl;

    std::cin.get();
    return 0;
}