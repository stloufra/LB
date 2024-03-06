//
// Created by stloufra on 10/30/23.
//

#ifndef INLETVELOCITY_H
#define INLETVELOCITY_H

#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct InletVelocity {
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void inlet(LBMDataPointer &Data, LBMConstantsPointer &Constants) {


        auto inlet_view = Data->meshBoundaryInlet.getView();

        auto rho_view = Data->rho.getView();
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();

        const auto Nvel = Constants->Nvel;

        MODELDATA MD;

        auto bb_inlet = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {
            Vertex vert = inlet_view[i.x()].vertex;
            Vector norm = inlet_view[i.x()].normal;
            Vector velc = inlet_view[i.x()].velocity;


            for (int vel = 0; vel < Nvel; vel++) {

                if (norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] > 0) {
                    df_view(vert.x, vert.y, vert.z, MD.c_rev[vel]) = df_post_view(vert.x, vert.y, vert.z, vel)
                                                                     - 6 * MD.weight[vel] *
                                                                       rho_view(vert.x, vert.y, vert.z) * (
                                                                               MD.c[vel][0] * velc.x() +
                                                                               MD.c[vel][1] * velc.y() +
                                                                               MD.c[vel][2] * velc.z()
                                                                       );

                }

            }

        };


        parallelFor<DeviceType>(0, Constants->inlet_num, bb_inlet);

    }

};

#endif
