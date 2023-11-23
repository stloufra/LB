//
// Created by stloufra on 10/30/23.
//

#ifndef OUTLETDENSITYEQUILIBRIUM_H
#define OUTLETDENSITYEQUILIBRIUM_H

#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct OutletDensityEquilibrium {
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void outlet(LBMDataPointer &Data, LBMConstantsPointer &Constants) {


        auto outlet_ConstView = Data->meshBoundaryOutlet.getConstView();
        auto u_ConstView = Data->u.getConstView();

        auto rho_ConstView = Data->rho.getConstView();
        auto df_post_ConstView = Data->df_post.getConstView();

        auto df_view = Data->df.getView();

        const auto Nvel = Constants->Nvel;

        MODELDATA MD;

        auto f_equilibrium_outlet = [=]
        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k,
        const int &vel,
        const RealType &density) mutable
        {
            RealType uc, u2;

            uc = MD.c[vel][0] * u_ConstView(i, j, k).x()
                 + MD.c[vel][1] * u_ConstView(i, j, k).y()
                 + MD.c[vel][2] * u_ConstView(i, j, k).z();

            u2 = u_ConstView(i, j, k).x() * u_ConstView(i, j, k).x()
                 + u_ConstView(i, j, k).y() * u_ConstView(i, j, k).y()
                 + u_ConstView(i, j, k).z() * u_ConstView(i, j, k).z();


            return MD.weight[vel] * density * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };

        auto bb_outlet = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {
            Vertex vert = outlet_ConstView[i.x()].vertex;
            Vector norm = outlet_ConstView[i.x()].normal;
            RealType density = outlet_ConstView[i.x()].density;

            for (int vel = 0; vel < Nvel; vel++) {

                if (norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] > 0) {
                    df_view(vert.x, vert.y, vert.z, MD.c_rev[vel]) = f_equilibrium_outlet(vert.x, vert.y, vert.z,
                                                                                          vel, density);

                }

            }

        };


        parallelFor<DeviceType>(0, Constants->outlet_num, bb_outlet);

    }

};

#endif
