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
        auto u_view = Data->u.getView();
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

        VectorType u_0 = Constants->VelocityInit;


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


        auto init_variables = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i ) mutable
        {
            if (mesh_view(i.x(), i.y(), i.z()) != 0) {
                rho_view(i.x(), i.y(), i.z()) = rho_0 / Cm * Cl * Cl * Cl;
                u_view(i.x(), i.y(), i.z()).x() = u_0.x() * Cu_inverse;
                u_view(i.x(), i.y(), i.z()).y() = u_0.y() * Cu_inverse;
                u_view(i.x(), i.y(), i.z()).z() = u_0.z() * Cu_inverse;

            } else {
                rho_view(i.x(), i.y(), i.z()) = 0.f;
                u_view(i.x(), i.y(), i.z()).x() = 0.f;
                u_view(i.x(), i.y(), i.z()).y() = 0.f;
                u_view(i.x(), i.y(), i.z()).z() = 0.f;
            }
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
        parallelFor<DeviceType>(begin1, end1, init_variables);


        parallelFor<DeviceType>(begin1, end1, init_df);




        auto inletVelocities = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {

            Vector velc = inlet_view[i.x()].velocity;
            Vertex vert = inlet_view[i.x()].vertex;

            u_view(vert.x, vert.y, vert.z).x() = velc.x();
            u_view(vert.x, vert.y, vert.z).y() = velc.y();
            u_view(vert.x, vert.y, vert.z).z() = velc.z();

        };


        parallelFor<DeviceType>(0, Constants->inlet_num, inletVelocities);
    }

};

#endif
