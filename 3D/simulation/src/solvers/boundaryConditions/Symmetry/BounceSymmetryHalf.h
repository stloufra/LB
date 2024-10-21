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
        auto mesh_view = Data->meshFluid.getView();

        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();

        const auto Nvel = Constants->Nvel;

        MODELDATA MD;

        VectorType norm_x(1,0,0);
        VectorType norm_y(0,1,0);
        VectorType norm_z(0,0,1);


        auto abs_cast = []
        __cuda_callable__ (const VectorType& v) {
            return TNL::Containers::StaticVector<3, unsigned int>(
                std::abs(v[0]), std::abs(v[1]), std::abs(v[2])
            );
        };


        auto bb_sym = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {


            Vertex vert = symmetry_view[i.x()].vertex;
            auto norm = symmetry_view[i.x()].normal;
            auto reg = symmetry_view[i.x()].regular;

            /*if (abs_cast(norm) == norm_x) {
                for(int vel = 0; vel < Nvel; vel++)
                {
                    if (norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0) { //pokud miri dovnitr
                        if(MD.c_sym_x[vel] != vel)
                        {
                            df_view(vert.x, vert.y, vert.z, vel) = df_post_view(vert.x, vert.y, vert.z, MD.c_sym_x[vel]);
                        }else
                        {
                            int dx, dy, dz;

                            dx = vert.x - MD.c[vel][0];
                            dy = vert.y - MD.c[vel][1];
                            dz = vert.z - MD.c[vel][2];

                            if (mesh_view(dx, dy, dz) == 0) {
                                df_view(vert.x, vert.y, vert.z, vel) = df_post_view(vert.x, vert.y, vert.z, MD.c_rev[vel]);
                            }

                        }

                    }
                }
            }
            else */if (abs_cast(norm) == norm_y) {
                for(int vel = 0; vel < Nvel; vel++)
                {
                    if (norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0) {   //pokud miri dovnitr
                        if(MD.c[vel][1] != 0)
                        {
                            df_view(vert.x, vert.y, vert.z, vel) = df_post_view(vert.x, vert.y, vert.z, MD.c_sym_y[vel]);
                        }

                        int dx, dy, dz;

                        dx = vert.x - MD.c[vel][0];
                        dy = vert.y - MD.c[vel][1];
                        dz = vert.z - MD.c[vel][2];

                        if (mesh_view(dx, dy, dz) == 0) {
                            df_view(vert.x, vert.y, vert.z, vel) = df_post_view(vert.x, vert.y, vert.z, MD.c_rev[vel]);
                        }

                    }
                }
            }
            /*else if (abs_cast(norm) == norm_z) {
                for(int vel = 0; vel < Nvel; vel++)
                {
                    if (norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0) { //pokud miri dovnitr
                        if(MD.c_sym_x[vel] != vel)
                        {
                            df_view(vert.x, vert.y, vert.z, vel) = df_post_view(vert.x, vert.y, vert.z, MD.c_sym_z[vel]);
                        }else
                        {
                            int dx, dy, dz;

                            dx = vert.x - MD.c[vel][0];
                            dy = vert.y - MD.c[vel][1];
                            dz = vert.z - MD.c[vel][2];

                            if (mesh_view(dx, dy, dz) == 0) {
                                df_view(vert.x, vert.y, vert.z, vel) = df_post_view(vert.x, vert.y, vert.z, MD.c_rev[vel]);
                            }

                        }
                    }
                }
            }*/
            else {
                printf("Fail. \n");
            }
        };


        parallelFor<DeviceType>(0, Constants->symmetry_num, bb_sym);
    }


};

#endif
