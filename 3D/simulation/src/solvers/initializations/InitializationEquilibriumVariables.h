//
// Created by stloufra on 10/30/23.
//

#ifndef ININITIALIZATIONEQUILIBRIUMVARIABLES_H
#define ININITIALIZATIONEQUILIBRIUMVARIABLES_H

#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct InitializationEquilibriumVariables {


    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void initialization(LBMDataPointer &Data, LBMConstantsPointer &Constants) {
        auto rho_view = Data->rho.getView();
        auto u_view = Data->u.getView();
        auto rho_host_view = Data->rho_out.getView();
        auto u_host_view = Data->u_out.getView();

        auto mesh_view = Data->meshFluid.getView();
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();
        auto p_view = Data->p.getView();

        auto Cl = Constants->Cl;
        auto Ct = Constants->Ct;
        auto Cm = Constants->Cm;
        auto Cu_inverse = Constants->Cu_inverse;
        auto Cu = Constants->Cu;
        const auto Nvel = Constants->Nvel;
        auto rho_0 = Constants->rho_fyz;


        MODELDATA MD;

        auto filename = Constants->InitFileName;


        std::ifstream inputFile("results/" + filename);

        if (!inputFile.is_open()) {
            std::cerr << "Failed to open the file for reading." << std::endl;
            return;
        }

        std::string line;


        // Skip lines until "LOOKUP_TABLE default"
        while (std::getline(inputFile, line)) {
            if (line.find("LOOKUP_TABLE default") != std::string::npos) {
                break;
            }
        }

        // Read data into "rho" array
        RealType rho_holder;
        for (int k = 0; k < Constants->dimZ_int; k++) {

            for (int j = 0; j < Constants->dimY_int; j++) {

                for (int i = 0; i < Constants->dimX_int; i++) {

                    if (!(inputFile >> rho_holder)) {
                        std::cerr << "Error reading data for rho." << std::endl;
                        inputFile.close();
                        return;
                    }
                    rho_host_view(i, j, k) = rho_holder;
                }
            }
        }


        // Skip lines until "VECTORS U double"
        while (std::getline(inputFile, line)) {
            if (line.find("VECTORS U double") != std::string::npos) {
                break;
            }
        }

        // Read data into "u" array of arrays
        RealType u_x, u_y, u_z;
        for (int k = 0; k < Constants->dimZ_int; k++) {

            for (int j = 0; j < Constants->dimY_int; j++) {

                for (int i = 0; i < Constants->dimX_int; i++) {

                    if (!(inputFile >> u_x >> u_y >> u_z)) {
                        std::cerr << "Error reading data for u." << std::endl;
                        inputFile.close();
                        return;
                    }

                    u_host_view(i, j, k).x() = u_x;
                    u_host_view(i, j, k).y() = u_y;
                    u_host_view(i, j, k).z() = u_z;
                }
            }
        }

        inputFile.close();

        u_view = u_host_view;
        rho_view = rho_host_view;


        auto f_equilibrium = [=]
        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k,
        const int &velo ) mutable
        {
            RealType uc, u2;

            uc = MD.c[velo][0] * u_view(i, j, k).x()
                 + MD.c[velo][1] * u_view(i, j, k).y()
                 + MD.c[velo][2] * u_view(i, j, k).z();

            u2 = u_view(i, j, k).x() * u_view(i, j, k).x() + u_view(i, j, k).y() * u_view(i, j, k).y() +
                 u_view(i, j, k).z() * u_view(i, j, k).z();


            return MD.weight[velo] * rho_view(i, j, k) * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };


        auto init_df = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i ) mutable
        {

            for (int vel = 0; vel < Nvel; vel++) {
                if (mesh_view(i.x(), i.y(), i.z()) != 0) {
                    df_view(i.x(), i.y(), i.z(), vel) = f_equilibrium(i.x(), i.y(), i.z(), vel);
                    df_post_view(i.x(), i.y(), i.z(), vel) = f_equilibrium(i.x(), i.y(), i.z(), vel);

                } else {
                    df_view(i.x(), i.y(), i.z(), vel) = 0.f;
                    df_post_view(i.x(), i.y(), i.z(), vel) = 0.f;
                }
            }
        };


        TNL::Containers::StaticArray<3, int> begin1{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end1{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin1, end1, init_df);
    }

};

#endif
