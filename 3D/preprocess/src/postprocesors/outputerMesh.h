
#include <fstream>
#include <iostream>
#include <string>
#include "../traits/LBMTraits.h"

#ifndef OUTPUTERMESH_H
#define OUTPUTERMESH_H

class outputerMesh {
    using RealType = LBMTraits::RealType;
    using VectorType = LBMTraits::VectorType;
    using DeviceType = LBMTraits::DeviceType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

public:

    static void MeshMatrixOut(LBMDataPointer &Data, LBMConstantsPointer &Constants, std::string filename) {


        auto mesh_view = Data->meshFluidHost.getView();

        std::ofstream outputFile("../simulation/results/" + filename + ".txt");

        if (!outputFile.is_open()) {
            std::cerr << "Failed to open the file for writing." << std::endl;
            return;
        }

        // Write the NDMatrix parameters to the file
        outputFile << "Dimensions: " << Constants->dimX_int << " "
                   << Constants->dimY_int << " "
                   << Constants->dimZ_int << " " << std::endl;
        outputFile << "Resolution Factor: " << Constants->resolution_factor << std::endl;
        outputFile << "Additional Factor: " << Constants->additional_factor << std::endl;
        outputFile << "BBmaxx: " << Constants->BBmaxx << std::endl;
        outputFile << "BBmaxy: " << Constants->BBmaxy << std::endl;
        outputFile << "BBmaxz: " << Constants->BBmaxz << std::endl;
        outputFile << "BBminx: " << Constants->BBminx << std::endl;
        outputFile << "BBminy: " << Constants->BBminy << std::endl;
        outputFile << "BBminz: " << Constants->BBminz << std::endl;

        // Write the NDMatrix data to the file
        for (int z = 0; z < Constants->dimZ_int; ++z) {
            for (int y = 0; y < Constants->dimY_int; ++y) {
                for (int x = 0; x < Constants->dimX_int; ++x) {
                    outputFile << mesh_view(x, y, z) << " ";
                }
                outputFile << std::endl; // Move to the next row
            }
            outputFile << std::endl; // Separate the 2D slices
        }

        outputFile.close();
    }


    static void MeshMatrixIn(LBMDataPointer &Data, LBMConstantsPointer &Constants, std::string filename) {
        // loads the mesh matrix from the file stored in /results

        std::ifstream inputFile("../simulation/results" +filename+".txt");
        auto mesh_view = Data->meshFluidHost.getView();

        if (!inputFile.is_open()) {
            std::cerr << "Failed to open the file for reading." << std::endl;
            return;
        }

        std::string line;


        while (std::getline(inputFile, line)) {
            if (line.find("Dimensions:") != std::string::npos) {
                int xDim, yDim, zDim;
                sscanf(line.c_str(), "Dimensions: %d %d %d", &xDim, &yDim, &zDim);
                Constants->dimX_int = xDim;
                Constants->dimY_int = yDim;
                Constants->dimZ_int = zDim;
                Data->meshFluidHost.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);
            } else if (line.find("Resolution Factor:") != std::string::npos) {
                RealType resolutionFactor;
                sscanf(line.c_str(), "Scale Factor: %lf", &resolutionFactor);
                Constants->resolution_factor = resolutionFactor;

            } else if (line.find("Additional Factor:") != std::string::npos) {
                RealType additionalFactor;
                sscanf(line.c_str(), "Addition Factor: %lf", &additionalFactor);
                Constants->additional_factor = additionalFactor;
            } else if (line.find("BBmaxx:") != std::string::npos) {
                sscanf(line.c_str(), "BBmaxx: %lf", &Constants->BBmaxx);
            } else if (line.find("BBmaxy:") != std::string::npos) {
                sscanf(line.c_str(), "BBmaxy: %lf", &Constants->BBmaxy);
            } else if (line.find("BBmaxz:") != std::string::npos) {
                sscanf(line.c_str(), "BBmaxz: %lf", &Constants->BBmaxz);
            } else if (line.find("BBminx:") != std::string::npos) {
                sscanf(line.c_str(), "BBminx: %lf", &Constants->BBminx);
            } else if (line.find("BBminy:") != std::string::npos) {
                sscanf(line.c_str(), "BBminy: %lf", &Constants->BBminy);
            } else if (line.find("BBminz:") != std::string::npos) {
                sscanf(line.c_str(), "BBminz: %lf", &Constants->BBminz);
            } else if (line.find_first_not_of(" \t\r\n") != std::string::npos) {
                // Read the NDMatrix data
                for (int z = 0; z < Constants->dimX_int; ++z) {
                    for (int y = 0; y < Constants->dimY_int; ++y) {
                        for (int x = 0; x < Constants->dimZ_int; ++x) {
                            inputFile >> mesh_view(x, y, z);
                        }
                    }
                }
            }
        }

        inputFile.close();
    }

};

#endif //OUTPUTERMATRIX_H
