//
// Created by stloufra on 10/30/23.
//

#ifndef COLLISIONSRT_H
#define COLLISIONSRT_H

#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct CollisionSRTLaminar
{

    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void collision(LBMDataPointer &Data, LBMConstantsPointer &Constants) {
        auto rho_view = Data->rho.getView();

        auto ux_view = Data->ux.getView();
        auto uy_view = Data->uy.getView();
        auto uz_view = Data->uz.getView();

        auto mesh_view = Data->meshFluid.getView();
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();


        auto omega = Constants->omega;


        MODELDATA MD;

        const auto Nvel = Constants->Nvel;

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

        auto coll = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i )mutable
        {

            if (mesh_view(i.x(), i.y(), i.z()) != 0) {
                auto ux = ux_view(i.x(), i.y(), i.z());
                auto uy = uy_view(i.x(), i.y(), i.z());
                auto uz = uz_view(i.x(), i.y(), i.z());
                auto rho= rho_view(i.x(), i.y(), i.z());

                for (int vel = 0; vel < Nvel; vel++) {

                    df_post_view(i.x(), i.y(), i.z(), vel) =
                            (1-omega)*df_view(i.x(), i.y(), i.z(), vel)
                            +f_equilibrium(ux,uy,uz,rho, vel) * omega;
                }
            }
        };


        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, coll);
    }

};

#endif
