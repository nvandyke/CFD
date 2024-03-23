
#include <iostream>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <string>

using namespace std;

extern "C" {
    void FVrun(string, string, int);
}


inline bool exists(const string& name) {
    ifstream f(name.c_str());
    return f.good();
}



int main(int argc, char* argv[]) {
    std::cout.precision(15);
    std::cout << std::fixed;
    string initialfile = "";
    string meshfile = "meshes\\bump0.gri";
    int order = 1;
    int cores = 1;

    while ((++argv)[0]) {
        if (argv[0][0] == '-') {
            switch (argv[0][1]) {


            case 'i':
                if (exists((++argv)[0])) {
                    initialfile = argv[0];
                    std::cout << "init file: " << initialfile << std::endl;
                }
                break;
            case 'm':
                if (exists((++argv)[0])) {
                    meshfile = argv[0];
                    std::cout << "mesh file: " << meshfile << std::endl;
                }
                break;
            case 'o':
                order = atoi((++argv)[0]);
                std::cout << "order: " << order << std::endl;
                break;
            case 'c':
                cores = atoi((++argv)[0]);
                std::cout << "cores: " << cores << std::endl;
                break;
            case 'h':
            default:
                std::cout << "options are:" << std::endl;
                std::cout << "\t-i: input file" << std::endl;
                std::cout << "\t-m: mesh file" << std::endl;
                std::cout << "\t-o: order" << std::endl;
                std::cout << "\t-c: num cores" << std::endl;
                std::cout << "\t-h: help menu" << std::endl;
                return 1;
                break;
            }
        }

    }

    omp_set_num_threads(cores);
    auto start_time = chrono::high_resolution_clock::now();
    FVrun(initialfile, meshfile, order);
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    std::cout << "took " << duration.count() << "s" << std::endl;

    std::cin.get();
    return 0;
}