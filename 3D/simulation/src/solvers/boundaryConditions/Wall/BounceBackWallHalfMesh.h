//
// Created by stloufra on 10/30/23.
//

#ifndef BOUNCEBACKWALLHALFMESH_H
#define BOUNCEBACKWALLHALFMESH_H

#include "../../../traits/LBMTraits.h"

template<typename MODELDATA>
struct BounceBackWallHalfMesh
{
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void slipVelocity(LBMDataPointer& Data, LBMConstantsPointer& Constants)
    {
     }

    static void bounceBackWall(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        auto df_view = Data->df.getView();

        auto df_post_ConstView = Data->df_post.getConstView();
        auto mesh_ConstView = Data->meshFluid.getConstView();

        const auto Nvel = Constants->Nvel;
        const auto dimX_int = Constants->dimX_int;
        const auto dimY_int = Constants->dimY_int;
        const auto dimZ_int = Constants->dimZ_int;

        MODELDATA MD;

        auto bb_wall = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i  ) mutable
        {

            if(mesh_ConstView(i.x(), i.y(), i.z()) == 2) {

                int dx, dy, dz;

                for (int vel = 0; vel < Nvel; vel++) {

                    dx = i.x() - MD.c[vel][0];
                    dy = i.y() - MD.c[vel][1];
                    dz = i.z() - MD.c[vel][2];

                    if (mesh_ConstView(dx, dy, dz) == 0) {
                        df_view(i.x(), i.y(), i.z(), vel) = df_post_ConstView(i.x(), i.y(), i.z(), MD.c_rev[vel]);
                    }

                }
            }

        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{dimX_int, dimY_int, dimZ_int};
        parallelFor<DeviceType>(begin, end, bb_wall);
    }

};

#endif
