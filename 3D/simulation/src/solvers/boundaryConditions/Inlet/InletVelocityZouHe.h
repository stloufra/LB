//
// Created by stloufra on 10/30/23.
//

#ifndef INLETVELOCITYZOUHEOLB_H
#define INLETVELOCITYZOUHEOLB_H

#include "../../../traits/LBMTraits.h"

template<typename MODELDATA>
struct InletVelocityZouHeOLB {

    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void inlet(LBMDataPointer &Data, LBMConstantsPointer &Constants) {


        auto inlet_view = Data->meshBoundaryInlet.getView();

        auto rho_view = Data->rho.getView();
        auto ux_view = Data->ux.getView();
        auto uy_view = Data->uy.getView();
        auto uz_view = Data->uz.getView();


        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();

        const auto Nvel = Constants->Nvel;

        MODELDATA MD;

        auto f_equilibrium_inlet = [=]
        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k,
        const int &vel,
        const VectorType &u) mutable
        {
            RealType uc, u2;

            uc = MD.c[vel][0] * u.x()
                 + MD.c[vel][1] * u.y()
                 + MD.c[vel][2] * u.z();

            u2 = u.x() * u.x()
                 + u.y() * u.y()
                 + u.z() * u.z();


            return MD.weight[vel] * rho_view(i, j, k)  * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };

        auto f_equilibrium = [=]
        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k,
        const int &vel ) mutable
        {
            RealType uc, u2;

            uc = MD.c[vel][0] * ux_view(i, j, k)
                 + MD.c[vel][1] * uy_view(i, j, k)
                 + MD.c[vel][2] * uz_view(i, j, k);

            u2 = ux_view(i, j, k) * ux_view(i, j, k)
                 + uy_view(i, j, k) * uy_view(i, j, k)
                 + uz_view(i, j, k) * uz_view(i, j, k);


            return MD.weight[vel] * rho_view(i, j, k) * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };



        auto bb_inlet = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {
            Vertex vert = inlet_view[i.x()].vertex;
            Vector norm = inlet_view[i.x()].normal;
            Vector u = inlet_view[i.x()].velocity;

            // f_i  = f'_i - f'_i^eq + f_i^eq

            for (int vel = 0; vel < Nvel; vel++) {

                if (norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] > 0) {

                    RealType fNeqRev = df_post_view(vert.x, vert.y, vert.z, vel) - f_equilibrium(vert.x, vert.y, vert.z, vel);

                    df_view(vert.x, vert.y, vert.z, MD.c_rev[vel]) = fNeqRev + f_equilibrium_inlet(vert.x, vert.y, vert.z, MD.c_rev[vel], u);

                }

            }
        };


        parallelFor<DeviceType>(0, Constants->inlet_num, bb_inlet);

    }

};

#endif
