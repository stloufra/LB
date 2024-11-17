//
// Created by stloufra on 10/30/23.
//

#ifndef PERIODICDELTAP_H
#define PERIODICDELTAP_H

#include "../../../traits/LBMTraits.h"

template<typename MODELDATA>
struct PeriodicDeltaP {
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void periodic(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        auto periodic_view = Data->meshBoundaryPeriodicDP.getView();
        auto mesh_view = Data->meshFluid.getView();

        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();

        auto ux_view = Data->ux.getView();
        auto uy_view = Data->uy.getView();
        auto uz_view = Data->uz.getView();
        auto rho_view = Data->rho.getView();

        const auto Nvel = Constants->Nvel;

        MODELDATA MD;

        VectorType norm_x(1,0,0);
        VectorType norm_y(0,1,0);
        VectorType norm_z(0,0,1);

        auto f_equilibrium_defined = [=]

        __cuda_callable__(

        const RealType &ux,

        const RealType &uy,

        const RealType &uz,

        const RealType &density,

        const int& vel
        )
        mutable
        {
            RealType uc, u2;

            uc = MD.c[vel][0] * ux + MD.c[vel][1] * uy + MD.c[vel][2] * uz;

            u2 = ux * ux + uy * uy + uz * uz;

            return MD.weight[vel] * density * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };


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
            auto DeltaRho = periodic_view[nod.x()].DeltaRho;

            int i = vert.x;
            int j = vert.y;
            int k = vert.z;




            if (abs_cast(norm) == norm_x) {

                auto rho = rho_view(perIdx,j,k);
                auto Drho = rho + DeltaRho;
                auto ux = ux_view(perIdx,j,k);
                auto uy = uy_view(perIdx,j,k);
                auto uz = uz_view(perIdx,j,k);

                for(int vel = 0; vel < Nvel; vel++)
                {
                    if(norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0)
                    {
                        auto df_post =df_post_view(perIdx,j,k,vel); // copy post collision
                        RealType feqN = f_equilibrium_defined( ux,uy,uz, rho, vel);
                        RealType feqDP = f_equilibrium_defined( ux,uy,uz, Drho, vel);

                        df_post = feqDP + (df_post - feqN);

                        int id;
                        int jd;
                        int kd;
                        id = i + MD.c[vel][0];
                        jd = j + MD.c[vel][1];
                        kd = k + MD.c[vel][2];

                        df_view(id, jd, kd, vel) = df_post;


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
                            else if ( norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] > 0 ){//miri ven -> periodicky

                                auto df_post =df_post_view(perIdx,j,k,MD.c_rev[vel]); // copy post collision
                                RealType feqN = f_equilibrium_defined( ux,uy,uz, rho, MD.c_rev[vel]);
                                RealType feqDP = f_equilibrium_defined( ux,uy,uz, rho, MD.c_rev[vel]);

                                df_post = feqDP + ( df_post - feqN);

                                df_view(i,j,k, vel) = df_post;
                            }
                        }
                    }

                }
            }
            else if (abs_cast(norm) == norm_y) {

                auto rho = rho_view(i,perIdx,k);
                auto Drho = rho + DeltaRho;
                auto ux = ux_view(i,perIdx,k);
                auto uy = uy_view(i,perIdx,k);
                auto uz = uz_view(i,perIdx,k);

                for(int vel = 0; vel < Nvel; vel++)
                {
                    if(norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0)
                    {
                        auto df_post =df_post_view(i,perIdx,k,vel); // copy post collision
                        RealType feqN = f_equilibrium_defined( ux,uy,uz, rho, vel);
                        RealType feqDP = f_equilibrium_defined( ux,uy,uz, Drho, vel);

                        df_post = feqDP + ( df_post - feqN );

                        int id;
                        int jd;
                        int kd;
                        id = i + MD.c[vel][0];
                        jd = j + MD.c[vel][1];
                        kd = k + MD.c[vel][2];

                        df_view(id, jd, kd, vel) = df_post;


                    }

                    if(!reg)
                    {
                        int dx, dy, dz;

                        dx = vert.x - MD.c[vel][0];
                        dz = vert.z - MD.c[vel][2];

                        if (mesh_view(dx, vert.y, dz) == 0)
                        {
                            if ( norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] <= 0 ){
                                df_view(i,j,k, vel) = df_post_view(i,j,k, MD.c_rev[vel]);
                            }
                            else if ( norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] > 0 ){
                                auto df_post =df_post_view(i,perIdx,k,MD.c_rev[vel]); // copy post collision
                                RealType feqN = f_equilibrium_defined( ux,uy,uz, rho, MD.c_rev[vel]);
                                RealType feqDP = f_equilibrium_defined( ux,uy,uz, rho, MD.c_rev[vel]);

                                df_post = feqDP + ( df_post - feqN );

                                df_view(i,j,k, vel) = df_post;
                            }
                        }
                    }
                }
            }
            else if (abs_cast(norm) == norm_z) {

                auto rho = rho_view(i,j,perIdx);
                auto Drho = rho + DeltaRho;
                auto ux = ux_view(i,j,perIdx);
                auto uy = uy_view(i,j,perIdx);
                auto uz = uz_view(i,j,perIdx);

                for(int vel = 0; vel < Nvel; vel++)
                {
                    if(norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0)
                    {
                        auto df_post =df_post_view(i,j,perIdx,vel); // copy post collision
                        RealType feqN = f_equilibrium_defined( ux,uy,uz, rho, vel);
                        RealType feqDP = f_equilibrium_defined( ux,uy,uz, Drho, vel);

                        df_post = feqDP + ( df_post - feqN );

                        int id;
                        int jd;
                        int kd;
                        id = i + MD.c[vel][0];
                        jd = j + MD.c[vel][1];
                        kd = k + MD.c[vel][2];

                        df_view(id, jd, kd, vel) = df_post;


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
                                auto df_post =df_post_view(i,j,perIdx,MD.c_rev[vel]); // copy post collision
                                RealType feqN = f_equilibrium_defined( ux,uy,uz, rho, MD.c_rev[vel]);
                                RealType feqDP = f_equilibrium_defined( ux,uy,uz, rho, MD.c_rev[vel]);

                                df_post = feqDP + ( df_post - feqN );

                                df_view(i,j,k, vel) = df_post;

                            }
                        }
                    }

                }
            }
            else {
                printf("Fail periodic DeltaP. \n");
            }
        };


        parallelFor<DeviceType>(0, Constants->periodicDP_num, bb_per);
    }

};

#endif
