//
// Created by stloufra on 10/30/23.
//

#ifndef PERIODICDELTAP_H
#define PERIODICDELTAP_H

#include "../../../traits/LBMTraits.h"

template <typename MODELDATA>
struct PeriodicDeltaP
{
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void periodic(LBMDataPointer& Data, LBMConstantsPointer& Constants)
    {
        auto periodic_view = Data->meshBoundaryPeriodicDP.getView();
        auto mesh_view = Data->meshFluid.getView();

        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();

        auto ux_view = Data->ux.getView();
        auto uy_view = Data->uy.getView();
        auto uz_view = Data->uz.getView();
        auto rho_view = Data->rho.getView();

        auto pmDview = Data->periodicMapDevice.getView();

        const auto Nvel = Constants->Nvel;

        MODELDATA MD;

        VectorType norm_x(1, 0, 0);
        VectorType norm_y(0, 1, 0);
        VectorType norm_z(0, 0, 1);

        auto f_equilibrium_defined = [=]__cuda_callable__(
                 const RealType& ux,
                 const RealType& uy,
                 const RealType& uz,
                 const RealType& density,
                 const int& vel)
        mutable
        {
            RealType uc, u2;

            uc = MD.c[vel][0] * ux + MD.c[vel][1] * uy + MD.c[vel][2] * uz;

            u2 = ux * ux + uy * uy + uz * uz;

            return MD.weight[vel] * density * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };


        auto abs_cast = []



        __cuda_callable__(const VectorType & v)
        {
            return TNL::Containers::StaticVector < 3, unsigned int > (
                std::abs(v[0]), std::abs(v[1]), std::abs(v[2])
            );
        };


        auto bb_per = [=]



        __cuda_callable__(



        const TNL::Containers::StaticArray<1, int>& nod
        )
        mutable
        {
            Vertex vert = periodic_view[nod.x()].vertex;
            auto norm = periodic_view[nod.x()].normal;
            auto reg = periodic_view[nod.x()].regular;
            auto perIdx = periodic_view[nod.x()].periodicIndex;
            auto DeltaRho = periodic_view[nod.x()].DeltaRho;

            int x = vert.x;
            int y = vert.y;
            int z = vert.z;

            //#define DEBUG
#ifdef DEBUG
            if(nod.x() == 4575)
            {
                printf("(x,y,z)=(%d,%d,%d) "
       "| norm - %f, %f, %f \n", vert.x, vert.y, vert.z, norm(0), norm(1), norm(2));
            }
#endif
#undef DEBUG

            if (abs_cast(norm) == norm_x)
            {
                norm = periodic_view[nod.x()].normal;

                auto rho = rho_view(perIdx, y, z);
                auto Drho = 1.f + DeltaRho;
                auto ux = ux_view(perIdx, y, z);
                auto uy = uy_view(perIdx, y, z);
                auto uz = uz_view(perIdx, y, z);

                for (int vel = 0; vel < Nvel; vel++)
                {
                    if (norm.x() * MD.c[vel][0] < 0)
                    {
                        auto df_post = df_post_view(perIdx, y, z, vel); // copy post collision
                        RealType feqN = f_equilibrium_defined(ux, uy, uz, rho, vel);
                        RealType feqDP = f_equilibrium_defined(ux, uy, uz, Drho, vel);

                        df_post = feqDP + (df_post - feqN);

                        int xd;
                        int yd;
                        int zd;
                        xd = x;
                        yd = y + MD.c[vel][1];
                        zd = z + MD.c[vel][2];

                        df_view(xd, yd, zd, vel) = df_post;
                    }

                    if (!reg) //BB Wall
                    {
                        int dy, dz;

                        dy = y + MD.c[vel][1];
                        dz = z + MD.c[vel][2];

                        auto xp = x - norm.x();

                        if (mesh_view(vert.x, dy, dz) == 0 && mesh_view(xp, y, z) != -3)
                        {
                            auto rhoNR = rho_view(x, y, z);
                            auto DrhoNR = 1.f + DeltaRho;
                            auto uxNR = ux_view(x, y, z);
                            auto uyNR = uy_view(x, y, z);
                            auto uzNR = uz_view(x, y, z);

                            auto df_post = df_post_view(x, y, z, vel);
                            RealType feqN = f_equilibrium_defined(uxNR, uyNR, uzNR, rhoNR, vel);
                            RealType feqDP = f_equilibrium_defined(uxNR, uyNR, uzNR, DrhoNR, vel);

                            df_post = feqDP + (df_post - feqN);


                            df_view(x, y, z, MD.c_rev[vel]) = df_post;
                        }
                    }
                }

                auto xp = x - norm.x();

                if(!reg){

                    if (mesh_view(xp, y, z) == -3)
                    {
                        //first y
                        VectorType yp(0, 1, 0);
                        VectorType ym(0, -1, 0);

                        if (y == 1)
                        {
                            for (int i = 0; i < pmDview.getSize(); i++)
                            {
                                VectorType normP = -pmDview[i].partner_normal;
                                int indexP = pmDview[i].partner_index;
                                int key = pmDview[i].key;

                                if (normP == ym && key == y)
                                {
                                    VectorType in = -(normP + norm);

                                    //#define DEBUG
#ifdef DEBUG
                                    if (nod.x() == 4575)
                                    {
                                        printf("(x,y,z)=(%d,%d,%d) "
                                               "| IN - %f, %f, %f "
                                               "| normP-(%f,%f,%f) "
                                               "| norm-(%f,%f,%f) "
                                               "| indexP - %d "
                                               "| key - %d \n", x, y, z, in(0), in(1), in(2), normP(0), normP(1),
                                               normP(2), norm(0), norm(1), norm(2), indexP, key);
                                    }
#endif
#undef DEBUG

                                    for (int vel = 0; vel < Nvel; vel++)
                                    {
                                        VectorType c(MD.c[vel][0], MD.c[vel][1], MD.c[vel][2]);

                                        if (in(0) == c(0) && in(1) == c(1))
                                        {


                                            auto rhoNR = rho_view(perIdx, indexP, z);
                                            auto DrhoNR = 1.f + DeltaRho;
                                            auto uxNR = ux_view(perIdx, indexP, z);
                                            auto uyNR = uy_view(perIdx, indexP, z);
                                            auto uzNR = uz_view(perIdx, indexP, z);

                                            auto df_post = df_post_view(perIdx, indexP, z, vel);

                                            RealType feqN = f_equilibrium_defined(uxNR, uyNR, uzNR, rhoNR, vel);
                                            RealType feqDP = f_equilibrium_defined(uxNR, uyNR, uzNR, DrhoNR, vel);

                                            df_post = feqDP + (df_post - feqN);

//#define DEBUG
#ifdef DEBUG
                                            if (nod.x() == 4575){
                                                printf(" IN - %f, %f, %f "
                                               "| c-(%f,%f,%f) "
                                               "|from - (%d,%d,%d)"
                                               "|to - (%d,%d,%d) \n", in(0), in(1), in(2), c(0), c(1), c(2),perIdx, indexP, z, x, y, z);
                                            }
#endif
#undef DEBUG


                                            df_view(x, y, z, vel) = df_post;
                                        }
                                    }
                                }
                            }
                        }
                        if (y == (Constants->dimY_int - 2))
                        {
                            for (int i = 0; i < pmDview.getSize(); i++)
                            {
                                auto normP = -pmDview[x].partner_normal;
                                auto indexP = pmDview[x].partner_index;
                                auto key = pmDview[x].key;

                                if (normP == yp && key == y)
                                {
                                    VectorType in = -1 * (normP + norm);

                                    for (int vel = 0; vel < Nvel; vel++)
                                    {
                                        VectorType c(MD.c[vel][0], MD.c[vel][1], MD.c[vel][2]);

                                        if (in(0) == c(0) && in(1) == c(1))
                                        {
                                            auto rhoNR = rho_view(perIdx, indexP, z);
                                            auto DrhoNR = 1.f + DeltaRho;
                                            auto uxNR = ux_view(perIdx, indexP, z);
                                            auto uyNR = uy_view(perIdx, indexP, z);
                                            auto uzNR = uz_view(perIdx, indexP, z);

                                            auto df_post = df_post_view(perIdx, indexP, z, vel);

                                            RealType feqN = f_equilibrium_defined(uxNR, uyNR, uzNR, rhoNR, vel);
                                            RealType feqDP = f_equilibrium_defined(uxNR, uyNR, uzNR, DrhoNR, vel);

                                            df_post = feqDP + (df_post - feqN);


                                            df_view(x, y, z, vel) = df_post;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /*else if (abs_cast(norm) == norm_y) {

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
            }*/
            else
            {
                printf("Fail periodic DeltaP. \n");
            }
        };


        parallelFor<DeviceType>(0, Constants->periodicDP_num, bb_per);
    }
};

#endif
