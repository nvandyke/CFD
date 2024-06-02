
#include <iostream>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <string>


extern "C" {
    void FVrun(std::string, std::string, int, double, double);
}


inline bool exists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}



int main(int argc, char* argv[]) {
    std::cout.precision(15);
    std::cout << std::fixed;
    std::string initialfile = "";
    std::string meshfile = "meshes\\bump0.gri";
    int order = 1;
    int cores = 1;
    double mach = 0.5;
    double aoa = 0.0;

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
            case 'M':
                mach = atof((++argv)[0]);
                std::cout << "Mach Number: " << mach << std::endl;
                break;
            case 'a':
                aoa = atof((++argv)[0]);
                std::cout << "Angle of Attack: " << aoa << std::endl;
                break;
            case 'h':
            default:
                std::cout << "options are:" << std::endl;
                std::cout << "\t-i: input file" << std::endl;
                std::cout << "\t-m: mesh file" << std::endl;
                std::cout << "\t-o: order" << std::endl;
                std::cout << "\t-M: Mach Number" << std::endl;
                std::cout << "\t-a: Angle of Attack" << std::endl;
                std::cout << "\t-c: num cores" << std::endl;
                std::cout << "\t-h: help menu" << std::endl;
                return 1;
                break;
            }
        }

    }

    omp_set_num_threads(cores);
    auto start_time = std::chrono::high_resolution_clock::now();
    FVrun(initialfile, meshfile, order, mach, aoa);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    std::cout << "took " << duration.count() << "s" << std::endl;

    std::cin.get();
    return 0;
}