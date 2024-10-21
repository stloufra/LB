//
// Created by stloufra on 10/16/23.
//

#include <fstream>
#include <iostream>
#include <string>
#include "../traits/LBMTraits.h"

#ifndef INC_24_4_GEO_TO_TNL_OUTPUTER_H
#define INC_24_4_GEO_TO_TNL_OUTPUTER_H

class outputerVTK
{
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

public:

    static void MeshVTK(LBMDataPointer &Data, LBMConstantsPointer &Constants, std::string filename) {
        std::cout << "\nWriting vtk file - " << filename << ".vtk\n";

        std::ofstream out_file( "../simulation/results/" + filename + ".vtk");

        out_file << "# vtk DataFile Version 2.0\n";
        out_file << "LBE " << filename << "\n";
        out_file << "ASCII\n";
        out_file << "DATASET STRUCTURED_POINTS\n";
        out_file << "DIMENSIONS " << Constants->dimX_int << " " << Constants->dimY_int << " " << Constants->dimZ_int
                 << "\n";
        out_file << "ASPECT_RATIO " << 1.f / Constants->resolution_factor << " " << 1.f / Constants->resolution_factor
                 << " "
                 << 1.f / Constants->resolution_factor << "\n";
        out_file << "ORIGIN " << Constants->BBminx - (Constants->additional_factor-0.5) / Constants->resolution_factor<< " "
                 << Constants->BBminy - (Constants->additional_factor-0.5) / Constants->resolution_factor << " "
                 << Constants->BBminz - (Constants->additional_factor-0.5) /  Constants->resolution_factor<< "\n";
        out_file << "POINT_DATA " << Constants->dimX_int * Constants->dimY_int * Constants->dimZ_int << "\n";

        out_file << "SCALARS " << "object " << "double 1\n";
        out_file << "LOOKUP_TABLE default\n";

        for (int k = 0; k < Constants->dimZ_int; k++) {
            for (int j = 0; j < Constants->dimY_int; j++) {
                for (int i = 0; i < Constants->dimX_int; i++) {
                    out_file << Data->meshFluidHost(i, j, k) << "\n";
                }
            }
        }

        std::cout << "Mesh written.\n";
    }




};

#endif //INC_24_4_GEO_TO_TNL_OUTPUTER_H
