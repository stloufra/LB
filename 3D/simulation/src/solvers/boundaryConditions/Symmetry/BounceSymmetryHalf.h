//
// Created by stloufra on 10/30/23.
//

#ifndef BOUNCESYMMETRYHALF_H
#define BOUNCESYMMETRYHALF_H

#include "../../../traits/LBMTraits.h"

template<typename MODELDATA>
struct BounceSymmetryHalf {
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void symmetry(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        auto symmetry_view = Data->meshBoundarySymmetry.getView();

        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();

        const auto Nvel = Constants->Nvel;

        auto mesh_view = Data->meshFluid.getView();

        MODELDATA MD;

        auto bb_sym = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {


            Vertex vert = symmetry_view[i.x()].vertex;
            auto norm = symmetry_view[i.x()].normal;

            if (norm == 1) {
                for (int vel = 0; vel < Nvel; vel++) {

                    int dx, dy, dz;

                    dx = vert.x - MD.c[vel][0];
                    dy = vert.y - MD.c[vel][1];
                    dz = vert.z - MD.c[vel][2];

                    if (mesh_view(dx, dy, dz) == 0) {

                        df_view(vert.x, vert.y, vert.z, vel) = df_post_view(vert.x, vert.y, vert.z, MD.sym_x[vel]);
                    }

                }
            } else if (norm == 1) {
                for (int vel = 0; vel < Nvel; vel++) {
                    int dx, dy, dz;

                    dx = vert.x - MD.c[vel][0];
                    dy = vert.y - MD.c[vel][1];
                    dz = vert.z - MD.c[vel][2];

                    if (mesh_view(dx, dy, dz) == 0) {

                        df_view(vert.x, vert.y, vert.z, vel) = df_post_view(vert.x, vert.y, vert.z, MD.sym_y[vel]);
                    }
                }
            } else if (norm == 2) {
                for (int vel = 0; vel < Nvel; vel++) {

                    int dx, dy, dz;

                    dx = vert.x - MD.c[vel][0];
                    dy = vert.y - MD.c[vel][1];
                    dz = vert.z - MD.c[vel][2];

                    if (mesh_view(dx, dy, dz) == 0) {

                        df_view(vert.x, vert.y, vert.z, vel) = df_post_view(vert.x, vert.y, vert.z, MD.sym_z[vel]);
                    }
                }
            } else {
                printf("Fail. \n");
            }
        };


        parallelFor<DeviceType>(0, Constants->symmetry_num, bb_sym);
    }


};

#endif
