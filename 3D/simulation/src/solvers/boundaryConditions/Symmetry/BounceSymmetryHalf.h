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

        VectorType norm_x(1, 0, 0);
        VectorType norm_y(0, 1, 0);
        VectorType norm_z(0, 0, 1);


        auto abs_cast = []
        __cuda_callable__(const VectorType & v)
        {
            return TNL::Containers::StaticVector < 3, unsigned int > (
                std::abs(v[0]), std::abs(v[1]), std::abs(v[2])
            );
        };


        auto bb_sym = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int>& i
        )
        mutable
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
            else */
            if (abs_cast(norm) == norm_y)
            {
                int x = vert.x;
                int y = vert.y;
                int z = vert.z;

                for (int vel = 0; vel < Nvel; vel++)
                {
                    int dx = x - MD.c[vel][0];
                    int dz = z - MD.c[vel][2];


                    if(reg)
                    {
                        if (norm.y() * MD.c[vel][1]< 0){
                            //pokud miri dovnitr a je regularni
#ifdef DEBUG
                            if(vel != 22 && vel != 20)
                            {
                                printf("Got here. My idx: {%d, %d, %d}. My vel %d [%d,%d,%d]. My sym vel %d. Point I want to cope from[%d,%d,%d] \n"
                                    , x, y, z, vel, MD.c[vel][0], MD.c[vel][1], MD.c[vel][2], MD.c_sym_y[vel], dx, y, dz);
                            }
#endif
                            if(vel == 20)
                            {
                                printf("Fucked up.");
                            }
                            df_view(x, y, z, vel) = df_post_view(dx, y, dz, MD.c_sym_y[vel]);
                        }
                    }

                    /*else if (!reg)
                    {
                        if (mesh_view(dx, y, dz) != 0 && mesh_view(dx, y, dz) != -2)
                        {
                            //not outlet not wall
                            if (norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0)
                            {
                                df_view(x, y, z, vel) = df_post_view(dx, y, dz, MD.c_sym_y[vel]);
                            }
                        }
                        else if (mesh_view(dx, y, dz) == -2)
                        {

                            if (norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0)
                            {
                                //outlet
                                df_view(x, y, z, vel) = df_post_view(dx, y, dz, MD.c_sym_y[vel]);
                                // a jeste tam na kopírovat
                                df_view(dx, y, dz, MD.c_sym_y[MD.c_rev[vel]]) = df_post_view(dx, y, dz, MD.c_rev[vel]);
                            }


                        }
                        else if (mesh_view(dx, y, dz) == 0)
                        {
                            //o stenu BB
                            df_view(x, y, z, vel) = df_post_view(x, y, z, MD.c_rev[vel]);
                        }
                    }*/
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
            else
            {
                printf("Fail. \n");
            }
        };


        parallelFor<DeviceType>(0, Constants->symmetry_num, bb_sym);
    }
};

#endif
