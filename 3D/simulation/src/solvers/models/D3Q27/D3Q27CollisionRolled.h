//
// Created by stloufra on 10/30/23.
//

#ifndef D3Q27COLLISIONROLLED_H
#define D3Q27COLLISIONROLLED_H

#include "../../../traits/LBMTraits.h"

struct D3Q27CollisionRolled
{

    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void collision(LBMDataPointer &Data, LBMConstantsPointer &Constants) {
        auto rho_view = Data->rho.getView();
        auto u_view = Data->u.getView();
        auto mesh_view = Data->meshFluid.getView();
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();


        auto omega = Constants->omega;


        D3Q27Data Mdata;

        const auto Nvel = Mdata.numberOfDiscreteVelocities;

        auto f_equilibrium = [=]
        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k,
        const int &vel ) mutable
        {
            RealType uc, u2;

            uc = Mdata.c[vel][0] * u_view(i, j, k).x()
                 + Mdata.c[vel][1] * u_view(i, j, k).y()
                 + Mdata.c[vel][2] * u_view(i, j, k).z();

            u2 = u_view(i, j, k).x() * u_view(i, j, k).x()
                 + u_view(i, j, k).y() * u_view(i, j, k).y()
                 + u_view(i, j, k).z() * u_view(i, j, k).z();


            return Mdata.weight[vel] * rho_view(i, j, k) * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };

        auto coll = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i )mutable
        {

            if (mesh_view(i.x(), i.y(), i.z()) != 0) {

                for (int vel = 0; vel < Nvel; vel++) {

                    df_post_view(i.x(), i.y(), i.z(), vel) =
                            df_view(i.x(), i.y(), i.z(), vel) -
                            (df_view(i.x(), i.y(), i.z(), vel) - f_equilibrium(i.x(), i.y(), i.z(), vel)) * omega;
                }
            }
            /*else {
                for (int vel = 0; vel < Nvel; vel++) {
                    df_post_view(i.x(), i.y(), i.z(), vel) = 0;
                }
            }*/
        };


        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, coll);
    }

};

#endif
