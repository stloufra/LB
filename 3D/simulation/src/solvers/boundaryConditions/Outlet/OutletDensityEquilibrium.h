//
// Created by stloufra on 10/30/23.
//

#ifndef OUTLETDENSITYEQUILIBRIUM_H
#define OUTLETDENSITYEQUILIBRIUM_H

#include "../../../traits/LBMTraits.h"

template<typename MODELDATA>
struct OutletDensityEquilibrium {
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void outletOmega(LBMDataPointer &Data, LBMConstantsPointer &Constants) {}

    static void outlet(LBMDataPointer &Data, LBMConstantsPointer &Constants) {


        auto outlet_view = Data->meshBoundaryOutlet.getView();
        auto ux_view = Data->ux.getView();
        auto uy_view = Data->uy.getView();
        auto uz_view = Data->uz.getView();

        auto rho_view = Data->rho.getView();
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();

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

            auto ux = ux_view(i, j, k);
            auto uy = uy_view(i, j, k);
            auto uz = uz_view(i, j, k);
            
            uc = MD.c[vel][0] * ux
                 + MD.c[vel][1] * uy
                 + MD.c[vel][2] * uz;

            u2 = ux * ux
                 + uy * uy
                 + uz * uz;


            return MD.weight[vel] * density * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };

        auto bb_outlet = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {
            Vertex vert = outlet_view[i.x()].vertex;
            Vector norm = outlet_view[i.x()].normal;
            RealType density = outlet_view[i.x()].density;

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
