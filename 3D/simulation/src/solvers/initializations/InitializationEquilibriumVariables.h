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

        auto ux_view = Data->ux.getView();
        auto uy_view = Data->uy.getView();
        auto uz_view = Data->uz.getView();

        auto rho_host_view = Data->rho_out.getView();

        auto ux_host_view = Data->u_x_out.getView();
        auto uy_host_view = Data->u_y_out.getView();
        auto uz_host_view = Data->u_z_out.getView();

        auto mesh_view = Data->meshFluid.getView();
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();
        auto p_view = Data->p.getView();
        auto inlet_view = Data->meshBoundaryInlet.getView();

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

                    ux_host_view(i, j, k) = u_x;
                    uy_host_view(i, j, k) = u_y;
                    uz_host_view(i, j, k) = u_z;
                }
            }
        }

        inputFile.close();

        ux_view = ux_host_view;
        uy_view = uy_host_view;
        uz_view = uz_host_view;


        rho_view = rho_host_view;


        auto f_equilibrium = [=]
        __cuda_callable__(
        RealType ux,
        RealType uy,
        RealType uz,
        RealType rho,
        const int &vel ) mutable
        {
            RealType uc, u2;

            uc = MD.c[vel][0] * ux
                 + MD.c[vel][1] * uy
                 + MD.c[vel][2] * uz;

            u2 = ux * ux
                 + uy * uy
                 + uz * uz;


            return MD.weight[vel] * rho * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };

        auto inletVelocities = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {

            Vector velc = inlet_view[i.x()].velocity;
            Vertex vert = inlet_view[i.x()].vertex;

            ux_view(vert.x, vert.y, vert.z) = velc.x();
            uy_view(vert.x, vert.y, vert.z) = velc.y();
            uz_view(vert.x, vert.y, vert.z) = velc.z();

        };

        parallelFor<DeviceType>(0, Constants->inlet_num, inletVelocities);


        auto init_df = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i ) mutable
        {
            auto ux = ux_view(i.x(), i.y(), i.z());
            auto uy = uy_view(i.x(), i.y(), i.z());
            auto uz = uz_view(i.x(), i.y(), i.z());
            auto rho= rho_view(i.x(), i.y(), i.z());

            for (int vel = 0; vel < Nvel; vel++) {
                if (mesh_view(i.x(), i.y(), i.z()) != 0) {
                    df_view(i.x(), i.y(), i.z(), vel) = f_equilibrium(ux,uy,uz,rho, vel);
                    df_post_view(i.x(), i.y(), i.z(), vel) = f_equilibrium(ux,uy,uz,rho, vel);

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
