//
// Created by stloufra on 10/30/23.
//

#ifndef ININITIALIZATIONEQUILIBRIUMCONSTVECTOR_H
#define ININITIALIZATIONEQUILIBRIUMCONSTVECTOR_H

#include "../models/D3Q27/D3Q27.h"
#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct InitializationEquilibriumConstVector {

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
        auto p_view = Data->p.getView();

        auto rhoTimeAvg_view = Data->rhoTimeAvg.getView();

        auto uxTimeAvg_view = Data->uxTimeAvg.getView();
        auto uyTimeAvg_view = Data->uyTimeAvg.getView();
        auto uzTimeAvg_view = Data->uzTimeAvg.getView();

        auto mesh_view = Data->meshFluid.getView();

        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();

        auto inlet_view = Data->meshBoundaryInlet.getView();

        auto Cl = Constants->Cl;
        auto Ct = Constants->Ct;
        auto Cm = Constants->Cm;
        auto Cu_inverse = Constants->Cu_inverse;
        auto Cu = Constants->Cu;

        const auto Nvel = Constants->Nvel;

        auto rho_0 = Constants->rho_fyz;

        MODELDATA MD;

        VectorType u_0 = Constants->VelocityInit;

        Constants->TimeAvgCounter = 0;
        Constants->timeAveraged = false;


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


        auto init_variables = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i ) mutable
        {
            if (mesh_view(i.x(), i.y(), i.z()) != 0) {
                rho_view(i.x(), i.y(), i.z()) = rho_0 / Cm * Cl * Cl * Cl;
                ux_view(i.x(), i.y(), i.z()) = u_0.x() * Cu_inverse;
                uy_view(i.x(), i.y(), i.z()) = u_0.y() * Cu_inverse;
                uz_view(i.x(), i.y(), i.z()) = u_0.z() * Cu_inverse;

            } else {
                rho_view(i.x(), i.y(), i.z()) = 0.f;
                ux_view(i.x(), i.y(), i.z()) = 0.f;
                uy_view(i.x(), i.y(), i.z()) = 0.f;
                uz_view(i.x(), i.y(), i.z()) = 0.f;
            }

            rhoTimeAvg_view(i.x(), i.y(), i.z()) = 0.f;
            uxTimeAvg_view(i.x(), i.y(), i.z()) = 0.f;
            uyTimeAvg_view(i.x(), i.y(), i.z()) = 0.f;
            uzTimeAvg_view(i.x(), i.y(), i.z()) = 0.f;
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
        parallelFor<DeviceType>(begin1, end1, init_variables);

        parallelFor<DeviceType>(0, Constants->inlet_num, inletVelocities);

        parallelFor<DeviceType>(begin1, end1, init_df);




   
    }

};

#endif
