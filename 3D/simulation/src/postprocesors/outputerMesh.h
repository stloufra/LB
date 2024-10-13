#pragma once
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

        std::ofstream outputFile("results/" + filename + ".txt");

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


    static void MeshMatrixIn(LBMDataPointer &Data, LBMConstantsPointer &Constants, std::string filename, bool verbose) {
        // loads the mesh matrix from the file stored in /results

        std::ifstream inputFile("results/" + filename + ".txt");

        if (!inputFile.is_open()) {
            std::cerr << "Failed to open the file for reading." << std::endl;
            return;
        }

        std::string line;

        if (std::is_same<RealType, double>::value) {
            while (std::getline(inputFile, line)) {
                if (line.find("Dimensions:") != std::string::npos) {
                    int xDim, yDim, zDim;
                    sscanf(line.c_str(), "Dimensions: %d %d %d", &xDim, &yDim, &zDim);
                    Constants->dimX_int = xDim;
                    Constants->dimY_int = yDim;
                    Constants->dimZ_int = zDim;

                } else if (line.find("Resolution Factor: ") != std::string::npos) {
                    RealType resolutionFactor;
                    sscanf(line.c_str(), "Resolution Factor: %lf", &resolutionFactor);
                    Constants->resolution_factor = resolutionFactor;

                } else if (line.find("Additional Factor: ") != std::string::npos) {
                    RealType additionalFactor;
                    sscanf(line.c_str(), "Additional Factor: %lf", &additionalFactor);
                    Constants->additional_factor = additionalFactor;

                } else if (line.find("BBmaxx:") != std::string::npos) {
                    RealType bbMaxx;
                    sscanf(line.c_str(), "BBmaxx: %lf", &bbMaxx);
                    Constants->BBmaxx = bbMaxx;

                } else if (line.find("BBmaxy:") != std::string::npos) {
                    RealType bbMaxy;
                    sscanf(line.c_str(), "BBmaxy: %lf", &bbMaxy);
                    Constants->BBmaxy = bbMaxy;

                } else if (line.find("BBmaxz:") != std::string::npos) {
                    RealType bbmaxz;
                    sscanf(line.c_str(), "BBmaxz: %lf", &bbmaxz);
                    Constants->BBmaxz = bbmaxz;

                } else if (line.find("BBminx:") != std::string::npos) {
                    RealType bbminx;
                    sscanf(line.c_str(), "BBminx: %lf", &bbminx);
                    Constants->BBminx = bbminx;

                } else if (line.find("BBminy:") != std::string::npos) {
                    RealType bbminy;
                    sscanf(line.c_str(), "BBminy: %lf", &bbminy);
                    Constants->BBminy = bbminy;

                } else if (line.find("BBminz:") != std::string::npos) {
                    RealType bbminz;
                    sscanf(line.c_str(), "BBminz: %lf", &bbminz);
                    Constants->BBminz = bbminz;
                    break;
                }

            }
        }

        else if (std::is_same<RealType, float>::value) {
            while (std::getline(inputFile, line)) {
                if (line.find("Dimensions:") != std::string::npos) {
                    int xDim, yDim, zDim;
                    sscanf(line.c_str(), "Dimensions: %d %d %d", &xDim, &yDim, &zDim);
                    Constants->dimX_int = xDim;
                    Constants->dimY_int = yDim;
                    Constants->dimZ_int = zDim;

                } else if (line.find("Resolution Factor: ") != std::string::npos) {
                    RealType resolutionFactor;
                    sscanf(line.c_str(), "Resolution Factor: %f", &resolutionFactor);
                    Constants->resolution_factor = resolutionFactor;

                } else if (line.find("Additional Factor: ") != std::string::npos) {
                    RealType additionalFactor;
                    sscanf(line.c_str(), "Additional Factor: %f", &additionalFactor);
                    Constants->additional_factor = additionalFactor;

                } else if (line.find("BBmaxx:") != std::string::npos) {
                    RealType bbMaxx;
                    sscanf(line.c_str(), "BBmaxx: %f", &bbMaxx);
                    Constants->BBmaxx = bbMaxx;

                } else if (line.find("BBmaxy:") != std::string::npos) {
                    RealType bbMaxy;
                    sscanf(line.c_str(), "BBmaxy: %f", &bbMaxy);
                    Constants->BBmaxy = bbMaxy;

                } else if (line.find("BBmaxz:") != std::string::npos) {
                    RealType bbmaxz;
                    sscanf(line.c_str(), "BBmaxz: %f", &bbmaxz);
                    Constants->BBmaxz = bbmaxz;

                } else if (line.find("BBminx:") != std::string::npos) {
                    RealType bbminx;
                    sscanf(line.c_str(), "BBminx: %f", &bbminx);
                    Constants->BBminx = bbminx;

                } else if (line.find("BBminy:") != std::string::npos) {
                    RealType bbminy;
                    sscanf(line.c_str(), "BBminy: %f", &bbminy);
                    Constants->BBminy = bbminy;

                } else if (line.find("BBminz:") != std::string::npos) {
                    RealType bbminz;
                    sscanf(line.c_str(), "BBminz: %f", &bbminz);
                    Constants->BBminz = bbminz;
                    break;
                }

            }
        }


        if (verbose) {
            std::cout << "Resolution factor: " << Constants->resolution_factor << std::endl;
            std::cout << "Additional factor: " << Constants->additional_factor << std::endl;
            std::cout << "Dimension X: " << Constants->dimX_int << std::endl;
            std::cout << "Dimension Y: " << Constants->dimY_int << std::endl;
            std::cout << "Dimension Z: " << Constants->dimZ_int << std::endl;
            std::cout << "BBminx: " << Constants->BBminx << std::endl;
            std::cout << "BBminy: " << Constants->BBminy << std::endl;
            std::cout << "BBminz: " << Constants->BBminz << std::endl;
            std::cout << "BBmaxx: " << Constants->BBmaxx << std::endl;
            std::cout << "BBmaxy: " << Constants->BBmaxy << std::endl;
            std::cout << "BBmaxz: " << Constants->BBmaxz << std::endl;
        }

        Data->meshFluidHost.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);

        auto mesh_view = Data->meshFluidHost.getView();

        int holder;
        for (int k = 0; k < Constants->dimZ_int; k++) {
            for (int j = 0; j < Constants->dimY_int; j++) {
                for (int i = 0; i < Constants->dimX_int; i++) {
                    inputFile >> holder;

                    mesh_view(i, j, k) = holder;
                }
            }
        }


        inputFile.close();
    }

};

#endif //OUTPUTERMATRIX_H
