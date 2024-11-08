//
// Created by stloufra on 10/30/23.
//

#ifndef PERIODIC_H
#define PERIODIC_H

#include "../../../traits/LBMTraits.h"

template<typename MODELDATA>
struct Periodic {
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void periodic(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        auto periodic_view = Data->meshBoundaryPeriodic.getView();
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


        auto bb_per = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &nod  ) mutable
        {


            Vertex vert = periodic_view[nod.x()].vertex;
            auto norm = periodic_view[nod.x()].normal;
            auto reg = periodic_view[nod.x()].regular;
            auto perIdx = periodic_view[nod.x()].periodicIndex;

            int i = vert.x;
            int j = vert.y;
            int k = vert.z;


            if (abs_cast(norm) == norm_x) {

                int perIdxC = perIdx - norm.x();

                for(int vel = 0; vel < Nvel; vel++)
                {
                    if(norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0)
                    {

                            df_view(i,j,k,vel)=df_view(perIdxC,j,k,vel);
                    }

                    if(!reg)
                    {
                        int dy, dz;

                        dy = vert.y + MD.c[vel][1];
                        dz = vert.z + MD.c[vel][2];

                        if (mesh_view(vert.x, dy, dz) == 0)
                        {
                            if ( norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] <= 0 ){ //miri dovnitr -> vlastni
                                df_view(i,j,k, vel) = df_post_view(i,j,k, MD.c_rev[vel]);
                            }
                            else if ( norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] > 0 ){ //miri ven -> periodicky
                                df_view(i,j,k, vel) = df_post_view(perIdx,j,k, MD.c_rev[vel]);
                            }
                        }
                    }

                }
            }
            else if (abs_cast(norm) == norm_y) {

                int perIdxC = perIdx - norm.y();

                for(int vel = 0; vel < Nvel; vel++)
                {
                    if(norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0)
                    {


                            df_view(i,j,k,vel)=df_view(i,perIdxC,k,vel);


                    }

                    if(!reg)
                    {
                        int dx, dy, dz;

                        dx = vert.x - MD.c[vel][0];
                        dy = vert.y - MD.c[vel][1];
                        dz = vert.z - MD.c[vel][2];

                        if (mesh_view(dx, vert.y, dz) == 0)
                        {
                            if ( norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] <= 0 ){
                                df_view(i,j,k, vel) = df_post_view(i,j,k, MD.c_rev[vel]);
                            }
                            else if ( norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] > 0 ){
                                df_view(i,j,k, vel) = df_post_view(i,perIdx,k, MD.c_rev[vel]);
                            }
                        }
                    }
                }
            }
            else if (abs_cast(norm) == norm_z) {

                int perIdxC = perIdx - norm.z();

                for(int vel = 0; vel < Nvel; vel++)
                {
                    if(norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0)
                    {

                            df_view(i,j,k,vel)=df_view(i,j,perIdxC,vel);

                    }

                    if(!reg)
                    {
                        int dx, dy, dz;

                        dx = vert.x + MD.c[vel][0];
                        dy = vert.y + MD.c[vel][1];
                        dz = vert.z + MD.c[vel][2];

                        if (mesh_view(dx, dy, vert.z) == 0)
                        {
                            if ( norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] <= 0 ){ //miri dovnitr -> vlastni
                                df_view(i,j,k, vel) = df_post_view(i,j,k, MD.c_rev[vel]);
                            }
                            else if ( norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] > 0 ){ //miri ven -> periodicky
                                df_view(i,j,k, vel) = df_post_view(i,j,perIdx, MD.c_rev[vel]);
                            }
                        }
                    }

                }
            }
            else {
                printf("Fail periodic. \n");
            }
        };


        parallelFor<DeviceType>(0, Constants->periodic_num, bb_per);
    }


};

#endif
