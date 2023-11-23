//
// Created by stloufra on 10/30/23.
//

#ifndef BOUNCEBACKWALLHALFVECTOR_H
#define BOUNCEBACKWALLHALFVECTOR_H

#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct BounceBackWallHalfVector
{
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void bounceBackWall(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        auto wall_ConstView = Data->meshBoundaryWall.getConstView();
        auto df_post_ConstView = Data->df_post.getConstView();
        auto mesh_ConstView = Data->meshFluid.getConstView();

        auto df_View = Data->df.getView();

        const auto Nvel = Constants ->Nvel;

        MODELDATA MD;

        auto bb_wall = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {


            Vertex vert = wall_ConstView[i.x()].vertex;

            int dx, dy, dz;

            for (int vel = 0; vel < Nvel; vel++) {

                dx = vert.x - MD.c[vel][0];
                dy = vert.y - MD.c[vel][1];
                dz = vert.z - MD.c[vel][2];

                if (mesh_ConstView(dx, dy, dz) == 0) {
                    df_View(vert.x, vert.y, vert.z, vel) = df_post_ConstView(vert.x, vert.y, vert.z, MD.c_rev[vel]);
                }

            }

        };

        parallelFor<DeviceType>(0, Constants->wall_num, bb_wall);

    }

};

#endif
