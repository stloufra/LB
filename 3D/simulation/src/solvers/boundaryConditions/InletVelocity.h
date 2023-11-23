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



        auto inlet_ConstView = Data->meshBoundaryInlet.getConstView();
        auto rho_ConstView = Data->rho.getConstView();
        auto df_post_ConstView = Data->df_post.getConstView();

        auto df_View = Data->df.getView();

        const auto Nvel = Constants->Nvel;

        MODELDATA MD;

        auto bb_inlet = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {
            Vertex vert = inlet_ConstView[i.x()].vertex;
            Vector norm = inlet_ConstView[i.x()].normal;
            Vector velc = inlet_ConstView[i.x()].velocity;


            for (int vel = 0; vel < Nvel; vel++) {

                if (norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] > 0) {
                    df_View(vert.x, vert.y, vert.z, MD.c_rev[vel]) = df_post_ConstView(vert.x, vert.y, vert.z, vel)
                                                                     - 6 * MD.weight[vel] *
                                                                       rho_ConstView(vert.x, vert.y, vert.z) * (
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
