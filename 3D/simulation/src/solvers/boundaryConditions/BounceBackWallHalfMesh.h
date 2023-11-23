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

        auto wall_view = Data->meshBoundaryWall;
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();
        const auto Nvel = Constants ->Nvel;
        auto mesh_view = Data->meshFluid.getView();

        MODELDATA MD;

        auto bb_wall = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {


            Vertex vert = wall_view[i.x()].vertex;

            int dx, dy, dz;

            for (int vel = 0; vel < Nvel; vel++) {

                dx = vert.x - MD.c[vel][0];
                dy = vert.y - MD.c[vel][1];
                dz = vert.z - MD.c[vel][2];

                if (mesh_view(dx, dy, dz) == 0) {
                    df_view(vert.x, vert.y, vert.z, vel) = df_post_view(vert.x, vert.y, vert.z, MD.c_rev[vel]);
                }

            }

        };

        parallelFor<DeviceType>(0, Constants->wall_num, bb_wall);

    }

};

#endif
