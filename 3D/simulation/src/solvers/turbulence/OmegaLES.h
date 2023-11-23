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
        auto u_view = Data->u.getView();
        auto omega_view = Data->omega.getView();

        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();

        auto mesh_view = Data->meshFluid.getView();


        const auto Nvel = Constants->Nvel;
        const auto cs4 = Constants->cs4;
        const auto tau = Constants->tau;
        const auto Cl = Constants -> Cl;
        const auto Ct = Constants -> Ct;
        const auto CLES = Constants -> CLES;


        MODELDATA MD;


        auto f_equilibrium = [=]
        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k,
        const int &velo ) mutable
        {
            RealType uc, u2;

            uc = MD.c[velo][0] * u_view(i, j, k).x()
                 + MD.c[velo][1] * u_view(i, j, k).y()
                 + MD.c[velo][2] * u_view(i, j, k).z();

            u2 = u_view(i, j, k).x() * u_view(i, j, k).x() + u_view(i, j, k).y() * u_view(i, j, k).y() +
                 u_view(i, j, k).z() * u_view(i, j, k).z();


            return MD.weight[velo] * rho_view(i, j, k) * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };


        auto init_omega = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i ) mutable
        {
            if (mesh_view(i.x(), i.y(), i.z()) != 0) {
                TensorType pi;
                RealType PI = 0.f;

                for (int alpha = 0; alpha < 3; alpha++) {
                    for (int beta = 0; beta < 3; beta++) {
                        for (int vel = 0; vel < Nvel; vel++) {

                            pi(alpha)(beta) += MD.c[vel][alpha] * MD.c[vel][beta] * (df_view(i.x(), i.y(), i.z(), vel) -
                                                                                     f_equilibrium(i.x(), i.y(), i.z(),
                                                                                                   vel));

                        }
                        PI += pi(alpha)(beta)*pi(alpha)(beta);
                    }
                }

                PI = sqrt(PI);

                RealType tauLES;

                tauLES = tau * 0.5f +
                         sqrt(tau * tau + 2.f * CLES * Cl * Cl / (rho_view(i.x(), i.y(), i.z()) * cs4 * Ct * Ct) * PI) *
                         0.5f;

                omega_view(i.x(), i.y(), i.z()) = 1.f / tauLES;
            }
        };


        TNL::Containers::StaticArray<3, int> begin1{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end1{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin1, end1, init_omega);
    }


};

#endif
