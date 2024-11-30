//
// Created by stloufra on 10/30/23.
//

#ifndef OMEGALES_H
#define OMEGALES_H

#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct OmegaLES {


    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using TensorType = LBMTraits::TensorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void omega(LBMDataPointer &Data, LBMConstantsPointer &Constants) {
        auto rho_view = Data->rho.getView();

        auto ux_view = Data->ux.getView();
        auto uy_view = Data->uy.getView();
        auto uz_view = Data->uz.getView();

        auto omega_view = Data->omega.getView();

        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();

        auto mesh_view = Data->meshFluid.getView();


        const auto Nvel = Constants->Nvel;
        const auto cs4 = Constants->cs4;
        const auto tau = Constants->tau;
        const auto Cl = Constants->Cl;
        const auto Ct = Constants->Ct;
        const auto CLES = Constants->CLES;
        const auto nu = Constants -> ny;


        MODELDATA MD;


        auto f_equilibrium = [=]
        __cuda_callable__(
        RealType ux,
        RealType uy,
        RealType uz,
        RealType rho,
        const int &velo ) mutable
        {
            RealType uc, u2;

            uc = MD.c[velo][0] * ux
                 + MD.c[velo][1] * uy
                 + MD.c[velo][2] * uz;

            u2 = ux * ux + uy * uy +
                 uz * uz;


            return MD.weight[velo] * rho * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };


        auto omegaFunc = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i ) mutable
        {
            if (mesh_view(i.x(), i.y(), i.z()) != 0) {

                auto u_x = ux_view(i.x(), i.y(), i.z());
                auto u_y = uy_view(i.x(), i.y(), i.z());
                auto u_z = uz_view(i.x(), i.y(), i.z());

                auto rho = rho_view(i.x(), i.y(), i.z());

                TensorType pi(0.f);


                for (int vel = 0; vel < Nvel; vel++) {

                    RealType feq = f_equilibrium(u_x, u_y, u_z, rho, vel);
                    RealType df = df_view(i.x(), i.y(), i.z(), vel);

                    pi(0)(0) += MD.c[vel][0] * MD.c[vel][0] * (df - feq);
                    pi(0)(1) += MD.c[vel][0] * MD.c[vel][1] * (df - feq);
                    pi(0)(2) += MD.c[vel][0] * MD.c[vel][2] * (df - feq);

                    pi(1)(0) += MD.c[vel][1] * MD.c[vel][0] * (df - feq);
                    pi(1)(1) += MD.c[vel][1] * MD.c[vel][1] * (df - feq);
                    pi(1)(2) += MD.c[vel][1] * MD.c[vel][2] * (df - feq);

                    pi(2)(0) += MD.c[vel][2] * MD.c[vel][0] * (df - feq);
                    pi(2)(1) += MD.c[vel][2] * MD.c[vel][1] * (df - feq);
                    pi(2)(2) += MD.c[vel][2] * MD.c[vel][2] * (df - feq);

                }
//#define DEBUG
#ifdef DEBUG

                if( i.x() == 100 && i.y() == 1 && i.z() == 40 )
                {
                    printf("Ux, Uy, Uz: %f, %f, %f \n", u_x, u_y, u_z);
                    printf("rho: %f \n", rho);

                    printf("PI:\n");

                    for( int i : {0,1,2})
                    {
                        for( int j : {0,1,2})
                        {
                            printf( "%f, ",pi(i)(j));
                        }
                        printf("\n");
                    }

                    printf("\n");

                    for (int vel = 0; vel < Nvel; vel++)
                    {
                        printf("F%d = %f \n", vel, df_view(i.x(), i.y(), i.z(), vel));
                    }

                    printf("\n");
                    printf("Ct - %f \n", Ct);
                    printf("Nu - %f \n", nu);

                }
#endif

                RealType PI = pi(0)(0) * pi(0)(0);
                PI += pi(0)(1) * pi(0)(1);
                PI += pi(0)(2) * pi(0)(2);

                PI += pi(1)(0) * pi(1)(0);
                PI += pi(1)(1) * pi(1)(1);
                PI += pi(1)(2) * pi(1)(2);

                PI += pi(2)(0) * pi(2)(0);
                PI += pi(2)(1) * pi(2)(1);
                PI += pi(2)(2) * pi(2)(2);

                RealType tauLES;

                RealType rho_inv = 1.0f / rho;
                RealType tau_sq = tau * tau;
                RealType CLES_term = 18.f * CLES * rho_inv;

                PI = sqrt(PI);

                // Compute tauLES with one square root
                tauLES = tau * 0.5f + 0.5f * sqrt(tau_sq + CLES_term * PI);

                //tauLES = tau * 0.5f +
                //         sqrt(tau * tau + 18.f * CLES / (rho) * sqrt(PI)) * 0.5f;

                omega_view(i.x(), i.y(), i.z()) = 1.f / tauLES;
#ifdef DEBUG
                if( i.x() == 100 && i.y() == 1 && i.z() == 40 )
                {
                    printf("Omega_tur = %f\n", omega_view(i.x(), i.y(), i.z()));
                }
#endif
            }
        };


        TNL::Containers::StaticArray<3, int> begin1{1, 1, 1};
        TNL::Containers::StaticArray<3, int> end1{Constants->dimX_int-1, Constants->dimY_int-1, Constants->dimZ_int-1};
        parallelFor<DeviceType>(begin1, end1, omegaFunc);
    }


};

#endif
