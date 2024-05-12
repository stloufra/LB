//
// Created by stloufra on 10/30/23.
//

#ifndef COLLISIONSRTUNROLL_H
#define COLLISIONSRTUNROLL_H

#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct CollisionSRTunroll {

    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void collision(LBMDataPointer &Data, LBMConstantsPointer &Constants) {
        auto rho_view = Data->rho.getView();
        auto u_view = Data->u.getView();
        auto mesh_view = Data->meshFluid.getView();
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();


        auto omega = Constants->omega;


        MODELDATA MD;

        auto coll = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i )mutable
        {

            if (mesh_view(i.x(), i.y(), i.z()) != 0) {

                df_post_view(i.x(), i.y(), i.z(), 0) = df_view(i.x(), i.y(), i.z(), 0) -
                                                       (df_view(i.x(), i.y(), i.z(), 0) -
                                                        0.296296f * rho_view(i.x(), i.y(), i.z()) *
                                                        (1.f
                                                         + 3.f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                  + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                  + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                         + 4.5f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                           * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                              + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                              + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                         - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                   u_view(i.x(), i.y(), i.z()).x()
                                                                   + u_view(i.x(), i.y(), i.z()).y() *
                                                                     u_view(i.x(), i.y(), i.z()).y()
                                                                   + u_view(i.x(), i.y(), i.z()).z() *
                                                                     u_view(i.x(), i.y(), i.z()).z())
                                                        )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 1) = df_view(i.x(), i.y(), i.z(), 1) -
                                                       (df_view(i.x(), i.y(), i.z(), 1) -
                                                        0.074074f * rho_view(i.x(), i.y(), i.z()) *
                                                        (1.f
                                                         + 3.f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                  + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                  + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                         + 4.5f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                           * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                              + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                              + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                         - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                   u_view(i.x(), i.y(), i.z()).x()
                                                                   + u_view(i.x(), i.y(), i.z()).y() *
                                                                     u_view(i.x(), i.y(), i.z()).y()
                                                                   + u_view(i.x(), i.y(), i.z()).z() *
                                                                     u_view(i.x(), i.y(), i.z()).z())
                                                        )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 2) = df_view(i.x(), i.y(), i.z(), 2) -
                                                       (df_view(i.x(), i.y(), i.z(), 2) -
                                                        0.074074f * rho_view(i.x(), i.y(), i.z()) *
                                                        (1.f
                                                         + 3.f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                  + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                  + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                         + 4.5f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                           * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                              + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                              + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                         - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                   u_view(i.x(), i.y(), i.z()).x()
                                                                   + u_view(i.x(), i.y(), i.z()).y() *
                                                                     u_view(i.x(), i.y(), i.z()).y()
                                                                   + u_view(i.x(), i.y(), i.z()).z() *
                                                                     u_view(i.x(), i.y(), i.z()).z())
                                                        )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 3) = df_view(i.x(), i.y(), i.z(), 3) -
                                                       (df_view(i.x(), i.y(), i.z(), 3) -
                                                        0.074074f * rho_view(i.x(), i.y(), i.z()) *
                                                        (1.f
                                                         + 3.f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                  + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                  + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                         + 4.5f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                           * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                              + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                              + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                         - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                   u_view(i.x(), i.y(), i.z()).x()
                                                                   + u_view(i.x(), i.y(), i.z()).y() *
                                                                     u_view(i.x(), i.y(), i.z()).y()
                                                                   + u_view(i.x(), i.y(), i.z()).z() *
                                                                     u_view(i.x(), i.y(), i.z()).z())
                                                        )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 4) = df_view(i.x(), i.y(), i.z(), 4) -
                                                       (df_view(i.x(), i.y(), i.z(), 4) -
                                                        0.074074f * rho_view(i.x(), i.y(), i.z()) *
                                                        (1.f
                                                         + 3.f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                  + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                  + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                         + 4.5f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                           * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                              + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                              + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                         - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                   u_view(i.x(), i.y(), i.z()).x()
                                                                   + u_view(i.x(), i.y(), i.z()).y() *
                                                                     u_view(i.x(), i.y(), i.z()).y()
                                                                   + u_view(i.x(), i.y(), i.z()).z() *
                                                                     u_view(i.x(), i.y(), i.z()).z())
                                                        )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 5) = df_view(i.x(), i.y(), i.z(), 5) -
                                                       (df_view(i.x(), i.y(), i.z(), 5) -
                                                        0.074074f * rho_view(i.x(), i.y(), i.z()) *
                                                        (1.f
                                                         + 3.f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                  + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                  + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                         + 4.5f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                           * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                              + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                              + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                         - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                   u_view(i.x(), i.y(), i.z()).x()
                                                                   + u_view(i.x(), i.y(), i.z()).y() *
                                                                     u_view(i.x(), i.y(), i.z()).y()
                                                                   + u_view(i.x(), i.y(), i.z()).z() *
                                                                     u_view(i.x(), i.y(), i.z()).z())
                                                        )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 6) = df_view(i.x(), i.y(), i.z(), 6) -
                                                       (df_view(i.x(), i.y(), i.z(), 6) -
                                                        0.074074f * rho_view(i.x(), i.y(), i.z()) *
                                                        (1.f
                                                         + 3.f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                  + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                  + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                         + 4.5f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                           * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                              + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                              + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                         - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                   u_view(i.x(), i.y(), i.z()).x()
                                                                   + u_view(i.x(), i.y(), i.z()).y() *
                                                                     u_view(i.x(), i.y(), i.z()).y()
                                                                   + u_view(i.x(), i.y(), i.z()).z() *
                                                                     u_view(i.x(), i.y(), i.z()).z())
                                                        )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 7) = df_view(i.x(), i.y(), i.z(), 7) -
                                                       (df_view(i.x(), i.y(), i.z(), 7) -
                                                        0.018519f * rho_view(i.x(), i.y(), i.z()) *
                                                        (1.f
                                                         + 3.f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                  + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                  + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                         + 4.5f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                           * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                              + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                              + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                         - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                   u_view(i.x(), i.y(), i.z()).x()
                                                                   + u_view(i.x(), i.y(), i.z()).y() *
                                                                     u_view(i.x(), i.y(), i.z()).y()
                                                                   + u_view(i.x(), i.y(), i.z()).z() *
                                                                     u_view(i.x(), i.y(), i.z()).z())
                                                        )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 8) = df_view(i.x(), i.y(), i.z(), 8) -
                                                       (df_view(i.x(), i.y(), i.z(), 8) -
                                                        0.018519f * rho_view(i.x(), i.y(), i.z()) *
                                                        (1.f
                                                         + 3.f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                  + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                  + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                         + 4.5f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                           * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                              + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                              + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                         - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                   u_view(i.x(), i.y(), i.z()).x()
                                                                   + u_view(i.x(), i.y(), i.z()).y() *
                                                                     u_view(i.x(), i.y(), i.z()).y()
                                                                   + u_view(i.x(), i.y(), i.z()).z() *
                                                                     u_view(i.x(), i.y(), i.z()).z())
                                                        )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 9) = df_view(i.x(), i.y(), i.z(), 9) -
                                                       (df_view(i.x(), i.y(), i.z(), 9) -
                                                        0.018519f * rho_view(i.x(), i.y(), i.z()) *
                                                        (1.f
                                                         + 3.f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                  + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                  + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                         + 4.5f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                           * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                              + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                              + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                         - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                   u_view(i.x(), i.y(), i.z()).x()
                                                                   + u_view(i.x(), i.y(), i.z()).y() *
                                                                     u_view(i.x(), i.y(), i.z()).y()
                                                                   + u_view(i.x(), i.y(), i.z()).z() *
                                                                     u_view(i.x(), i.y(), i.z()).z())
                                                        )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 10) = df_view(i.x(), i.y(), i.z(), 10) -
                                                        (df_view(i.x(), i.y(), i.z(), 10) -
                                                         0.018519f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                               + 0 * u_view(i.x(), i.y(), i.z()).y()
                                                               + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 11) = df_view(i.x(), i.y(), i.z(), 11) -
                                                        (df_view(i.x(), i.y(), i.z(), 11) -
                                                         0.018519f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                               + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                               + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 12) = df_view(i.x(), i.y(), i.z(), 12) -
                                                        (df_view(i.x(), i.y(), i.z(), 12) -
                                                         0.018519f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                               + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                               + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 13) = df_view(i.x(), i.y(), i.z(), 13) -
                                                        (df_view(i.x(), i.y(), i.z(), 13) -
                                                         0.018519f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                               + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                               + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 14) = df_view(i.x(), i.y(), i.z(), 14) -
                                                        (df_view(i.x(), i.y(), i.z(), 14) -
                                                         0.018519f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                               + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                               + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 15) = df_view(i.x(), i.y(), i.z(), 15) -
                                                        (df_view(i.x(), i.y(), i.z(), 15) -
                                                         0.018519f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                               + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                               + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 16) = df_view(i.x(), i.y(), i.z(), 16) -
                                                        (df_view(i.x(), i.y(), i.z(), 16) -
                                                         0.018519f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                               + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                               + 0 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 17) = df_view(i.x(), i.y(), i.z(), 17) -
                                                        (df_view(i.x(), i.y(), i.z(), 17) -
                                                         0.018519f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                               + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                               + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 18) = df_view(i.x(), i.y(), i.z(), 18) -
                                                        (df_view(i.x(), i.y(), i.z(), 18) -
                                                         0.018519f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (0 * u_view(i.x(), i.y(), i.z()).x()
                                                               + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                               + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 19) = df_view(i.x(), i.y(), i.z(), 19) -
                                                        (df_view(i.x(), i.y(), i.z(), 19) -
                                                         0.004630f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                               + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                               + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 20) = df_view(i.x(), i.y(), i.z(), 20) -
                                                        (df_view(i.x(), i.y(), i.z(), 20) -
                                                         0.004630f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                               + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                               + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 21) = df_view(i.x(), i.y(), i.z(), 21) -
                                                        (df_view(i.x(), i.y(), i.z(), 21) -
                                                         0.004630f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                               + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                               + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 22) = df_view(i.x(), i.y(), i.z(), 22) -
                                                        (df_view(i.x(), i.y(), i.z(), 22) -
                                                         0.004630f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                               + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                               + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 23) = df_view(i.x(), i.y(), i.z(), 23) -
                                                        (df_view(i.x(), i.y(), i.z(), 23) -
                                                         0.004630f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                               + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                               + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 24) = df_view(i.x(), i.y(), i.z(), 24) -
                                                        (df_view(i.x(), i.y(), i.z(), 24) -
                                                         0.004630f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                               + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                               + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 25) = df_view(i.x(), i.y(), i.z(), 25) -
                                                        (df_view(i.x(), i.y(), i.z(), 25) -
                                                         0.004630f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (-1 * u_view(i.x(), i.y(), i.z()).x()
                                                               + -1 * u_view(i.x(), i.y(), i.z()).y()
                                                               + -1 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;
                df_post_view(i.x(), i.y(), i.z(), 26) = df_view(i.x(), i.y(), i.z(), 26) -
                                                        (df_view(i.x(), i.y(), i.z(), 26) -
                                                         0.004630f * rho_view(i.x(), i.y(), i.z()) *
                                                         (1.f
                                                          + 3.f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                   + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                          + 4.5f * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                                    + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                                    + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                            * (1 * u_view(i.x(), i.y(), i.z()).x()
                                                               + 1 * u_view(i.x(), i.y(), i.z()).y()
                                                               + 1 * u_view(i.x(), i.y(), i.z()).z())
                                                          - 1.5f * (u_view(i.x(), i.y(), i.z()).x() *
                                                                    u_view(i.x(), i.y(), i.z()).x()
                                                                    + u_view(i.x(), i.y(), i.z()).y() *
                                                                      u_view(i.x(), i.y(), i.z()).y()
                                                                    + u_view(i.x(), i.y(), i.z()).z() *
                                                                      u_view(i.x(), i.y(), i.z()).z())
                                                         )) * omega;


            } else {

                df_post_view(i.x(), i.y(), i.z(), 0) = 0;
                df_post_view(i.x(), i.y(), i.z(), 1) = 0;
                df_post_view(i.x(), i.y(), i.z(), 2) = 0;
                df_post_view(i.x(), i.y(), i.z(), 3) = 0;
                df_post_view(i.x(), i.y(), i.z(), 4) = 0;
                df_post_view(i.x(), i.y(), i.z(), 5) = 0;
                df_post_view(i.x(), i.y(), i.z(), 6) = 0;
                df_post_view(i.x(), i.y(), i.z(), 7) = 0;
                df_post_view(i.x(), i.y(), i.z(), 8) = 0;
                df_post_view(i.x(), i.y(), i.z(), 9) = 0;
                df_post_view(i.x(), i.y(), i.z(), 10) = 0;
                df_post_view(i.x(), i.y(), i.z(), 11) = 0;
                df_post_view(i.x(), i.y(), i.z(), 12) = 0;
                df_post_view(i.x(), i.y(), i.z(), 13) = 0;
                df_post_view(i.x(), i.y(), i.z(), 14) = 0;
                df_post_view(i.x(), i.y(), i.z(), 15) = 0;
                df_post_view(i.x(), i.y(), i.z(), 16) = 0;
                df_post_view(i.x(), i.y(), i.z(), 17) = 0;
                df_post_view(i.x(), i.y(), i.z(), 18) = 0;
                df_post_view(i.x(), i.y(), i.z(), 19) = 0;
                df_post_view(i.x(), i.y(), i.z(), 20) = 0;
                df_post_view(i.x(), i.y(), i.z(), 21) = 0;
                df_post_view(i.x(), i.y(), i.z(), 22) = 0;
                df_post_view(i.x(), i.y(), i.z(), 23) = 0;
                df_post_view(i.x(), i.y(), i.z(), 24) = 0;
                df_post_view(i.x(), i.y(), i.z(), 25) = 0;
                df_post_view(i.x(), i.y(), i.z(), 26) = 0;


            }
        };


        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, coll);
    }

};

#endif
