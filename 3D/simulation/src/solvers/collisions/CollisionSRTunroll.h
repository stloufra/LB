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
        
      auto ux_view = Data->ux.getView();
      auto uy_view = Data->uy.getView();
      auto uz_view = Data->uz.getView();
        
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

              auto ux = ux_view(i.x(), i.y(), i.z());
              auto uy = uy_view(i.x(), i.y(), i.z());
              auto uz = uz_view(i.x(), i.y(), i.z());
              
              auto rho = rho_view(i.x(), i.y(), i.z());

                df_post_view(i.x(), i.y(), i.z(), 0) = (1-omega)*df_view(i.x(), i.y(), i.z(), 0)
                                                       +
                                                        0.296296f * rho *
                                                        (1.f
                                                         + 3.f * (0 * ux
                                                                  + 0 * uy
                                                                  + 0 * uz)
                                                         + 4.5f * (0 * ux
                                                                   + 0 * uy
                                                                   + 0 * uz)
                                                           * (0 * ux
                                                              + 0 * uy
                                                              + 0 * uz)
                                                         - 1.5f * (ux *
                                                                   ux
                                                                   + uy *
                                                                     uy
                                                                   + uz *
                                                                     uz)
                                                        ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 1) = (1-omega)*df_view(i.x(), i.y(), i.z(), 1)
                                                       +
                                                        0.074074f * rho *
                                                        (1.f
                                                         + 3.f * (1 * ux
                                                                  + 0 * uy
                                                                  + 0 * uz)
                                                         + 4.5f * (1 * ux
                                                                   + 0 * uy
                                                                   + 0 * uz)
                                                           * (1 * ux
                                                              + 0 * uy
                                                              + 0 * uz)
                                                         - 1.5f * (ux *
                                                                   ux
                                                                   + uy *
                                                                     uy
                                                                   + uz *
                                                                     uz)
                                                        ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 2) = (1-omega)*df_view(i.x(), i.y(), i.z(), 2)
                                                       +
                                                        0.074074f * rho *
                                                        (1.f
                                                         + 3.f * (-1 * ux
                                                                  + 0 * uy
                                                                  + 0 * uz)
                                                         + 4.5f * (-1 * ux
                                                                   + 0 * uy
                                                                   + 0 * uz)
                                                           * (-1 * ux
                                                              + 0 * uy
                                                              + 0 * uz)
                                                         - 1.5f * (ux *
                                                                   ux
                                                                   + uy *
                                                                     uy
                                                                   + uz *
                                                                     uz)
                                                        ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 3) = (1-omega)*df_view(i.x(), i.y(), i.z(), 3)
                                                       +
                                                        0.074074f * rho *
                                                        (1.f
                                                         + 3.f * (0 * ux
                                                                  + 0 * uy
                                                                  + -1 * uz)
                                                         + 4.5f * (0 * ux
                                                                   + 0 * uy
                                                                   + -1 * uz)
                                                           * (0 * ux
                                                              + 0 * uy
                                                              + -1 * uz)
                                                         - 1.5f * (ux *
                                                                   ux
                                                                   + uy *
                                                                     uy
                                                                   + uz *
                                                                     uz)
                                                        ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 4) = (1-omega)*df_view(i.x(), i.y(), i.z(), 4)
                                                       +
                                                        0.074074f * rho *
                                                        (1.f
                                                         + 3.f * (0 * ux
                                                                  + 0 * uy
                                                                  + 1 * uz)
                                                         + 4.5f * (0 * ux
                                                                   + 0 * uy
                                                                   + 1 * uz)
                                                           * (0 * ux
                                                              + 0 * uy
                                                              + 1 * uz)
                                                         - 1.5f * (ux *
                                                                   ux
                                                                   + uy *
                                                                     uy
                                                                   + uz *
                                                                     uz)
                                                        ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 5) = (1-omega)*df_view(i.x(), i.y(), i.z(), 5)
                                                       +
                                                        0.074074f * rho *
                                                        (1.f
                                                         + 3.f * (0 * ux
                                                                  + -1 * uy
                                                                  + 0 * uz)
                                                         + 4.5f * (0 * ux
                                                                   + -1 * uy
                                                                   + 0 * uz)
                                                           * (0 * ux
                                                              + -1 * uy
                                                              + 0 * uz)
                                                         - 1.5f * (ux *
                                                                   ux
                                                                   + uy *
                                                                     uy
                                                                   + uz *
                                                                     uz)
                                                        ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 6) = (1-omega)*df_view(i.x(), i.y(), i.z(), 6)
                                                       +
                                                        0.074074f * rho *
                                                        (1.f
                                                         + 3.f * (0 * ux
                                                                  + 1 * uy
                                                                  + 0 * uz)
                                                         + 4.5f * (0 * ux
                                                                   + 1 * uy
                                                                   + 0 * uz)
                                                           * (0 * ux
                                                              + 1 * uy
                                                              + 0 * uz)
                                                         - 1.5f * (ux *
                                                                   ux
                                                                   + uy *
                                                                     uy
                                                                   + uz *
                                                                     uz)
                                                        ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 7) = (1-omega)*df_view(i.x(), i.y(), i.z(), 7)
                                                       +
                                                        0.018519f * rho *
                                                        (1.f
                                                         + 3.f * (1 * ux
                                                                  + 0 * uy
                                                                  + -1 * uz)
                                                         + 4.5f * (1 * ux
                                                                   + 0 * uy
                                                                   + -1 * uz)
                                                           * (1 * ux
                                                              + 0 * uy
                                                              + -1 * uz)
                                                         - 1.5f * (ux *
                                                                   ux
                                                                   + uy *
                                                                     uy
                                                                   + uz *
                                                                     uz)
                                                        ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 8) = (1-omega)*df_view(i.x(), i.y(), i.z(), 8)
                                                       +
                                                        0.018519f * rho *
                                                        (1.f
                                                         + 3.f * (-1 * ux
                                                                  + 0 * uy
                                                                  + 1 * uz)
                                                         + 4.5f * (-1 * ux
                                                                   + 0 * uy
                                                                   + 1 * uz)
                                                           * (-1 * ux
                                                              + 0 * uy
                                                              + 1 * uz)
                                                         - 1.5f * (ux *
                                                                   ux
                                                                   + uy *
                                                                     uy
                                                                   + uz *
                                                                     uz)
                                                        ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 9) = (1-omega)*df_view(i.x(), i.y(), i.z(), 9)
                                                       +
                                                        0.018519f * rho *
                                                        (1.f
                                                         + 3.f * (1 * ux
                                                                  + 0 * uy
                                                                  + 1 * uz)
                                                         + 4.5f * (1 * ux
                                                                   + 0 * uy
                                                                   + 1 * uz)
                                                           * (1 * ux
                                                              + 0 * uy
                                                              + 1 * uz)
                                                         - 1.5f * (ux *
                                                                   ux
                                                                   + uy *
                                                                     uy
                                                                   + uz *
                                                                     uz)
                                                        ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 10) = (1-omega)*df_view(i.x(), i.y(), i.z(), 10)
                                                        +
                                                         0.018519f * rho *
                                                         (1.f
                                                          + 3.f * (-1 * ux
                                                                   + 0 * uy
                                                                   + -1 * uz)
                                                          + 4.5f * (-1 * ux
                                                                    + 0 * uy
                                                                    + -1 * uz)
                                                            * (-1 * ux
                                                               + 0 * uy
                                                               + -1 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 11) = (1-omega)*df_view(i.x(), i.y(), i.z(), 11)
                                                        +
                                                         0.018519f * rho *
                                                         (1.f
                                                          + 3.f * (-1 * ux
                                                                   + -1 * uy
                                                                   + 0 * uz)
                                                          + 4.5f * (-1 * ux
                                                                    + -1 * uy
                                                                    + 0 * uz)
                                                            * (-1 * ux
                                                               + -1 * uy
                                                               + 0 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 12) = (1-omega)*df_view(i.x(), i.y(), i.z(), 12)
                                                        +
                                                         0.018519f * rho *
                                                         (1.f
                                                          + 3.f * (1 * ux
                                                                   + 1 * uy
                                                                   + 0 * uz)
                                                          + 4.5f * (1 * ux
                                                                    + 1 * uy
                                                                    + 0 * uz)
                                                            * (1 * ux
                                                               + 1 * uy
                                                               + 0 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 13) = (1-omega)*df_view(i.x(), i.y(), i.z(), 13)
                                                        +
                                                         0.018519f * rho *
                                                         (1.f
                                                          + 3.f * (0 * ux
                                                                   + 1 * uy
                                                                   + -1 * uz)
                                                          + 4.5f * (0 * ux
                                                                    + 1 * uy
                                                                    + -1 * uz)
                                                            * (0 * ux
                                                               + 1 * uy
                                                               + -1 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 14) = (1-omega)*df_view(i.x(), i.y(), i.z(), 14) -
                                                        +
                                                         0.018519f * rho *
                                                         (1.f
                                                          + 3.f * (0 * ux
                                                                   + -1 * uy
                                                                   + 1 * uz)
                                                          + 4.5f * (0 * ux
                                                                    + -1 * uy
                                                                    + 1 * uz)
                                                            * (0 * ux
                                                               + -1 * uy
                                                               + 1 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 15) = (1-omega)*df_view(i.x(), i.y(), i.z(), 15)
                                                        +
                                                         0.018519f * rho *
                                                         (1.f
                                                          + 3.f * (-1 * ux
                                                                   + 1 * uy
                                                                   + 0 * uz)
                                                          + 4.5f * (-1 * ux
                                                                    + 1 * uy
                                                                    + 0 * uz)
                                                            * (-1 * ux
                                                               + 1 * uy
                                                               + 0 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 16) = (1-omega)*df_view(i.x(), i.y(), i.z(), 16)
                                                        +
                                                         0.018519f * rho *
                                                         (1.f
                                                          + 3.f * (1 * ux
                                                                   + -1 * uy
                                                                   + 0 * uz)
                                                          + 4.5f * (1 * ux
                                                                    + -1 * uy
                                                                    + 0 * uz)
                                                            * (1 * ux
                                                               + -1 * uy
                                                               + 0 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 17) = (1-omega)*df_view(i.x(), i.y(), i.z(), 17)
                                                        +
                                                         0.018519f * rho *
                                                         (1.f
                                                          + 3.f * (0 * ux
                                                                   + 1 * uy
                                                                   + 1 * uz)
                                                          + 4.5f * (0 * ux
                                                                    + 1 * uy
                                                                    + 1 * uz)
                                                            * (0 * ux
                                                               + 1 * uy
                                                               + 1 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 18) = (1-omega)*df_view(i.x(), i.y(), i.z(), 18)
                                                        +
                                                         0.018519f * rho *
                                                         (1.f
                                                          + 3.f * (0 * ux
                                                                   + -1 * uy
                                                                   + -1 * uz)
                                                          + 4.5f * (0 * ux
                                                                    + -1 * uy
                                                                    + -1 * uz)
                                                            * (0 * ux
                                                               + -1 * uy
                                                               + -1 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 19) = (1-omega)*df_view(i.x(), i.y(), i.z(), 19)
                                                        +
                                                         0.004630f * rho *
                                                         (1.f
                                                          + 3.f * (-1 * ux
                                                                   + 1 * uy
                                                                   + -1 * uz)
                                                          + 4.5f * (-1 * ux
                                                                    + 1 * uy
                                                                    + -1 * uz)
                                                            * (-1 * ux
                                                               + 1 * uy
                                                               + -1 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 20) = (1-omega)*df_view(i.x(), i.y(), i.z(), 20)
                                                        +
                                                         0.004630f * rho *
                                                         (1.f
                                                          + 3.f * (1 * ux
                                                                   + -1 * uy
                                                                   + 1 * uz)
                                                          + 4.5f * (1 * ux
                                                                    + -1 * uy
                                                                    + 1 * uz)
                                                            * (1 * ux
                                                               + -1 * uy
                                                               + 1 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 21) = (1-omega)*df_view(i.x(), i.y(), i.z(), 21)
                                                        +
                                                         0.004630f * rho *
                                                         (1.f
                                                          + 3.f * (-1 * ux
                                                                   + -1 * uy
                                                                   + 1 * uz)
                                                          + 4.5f * (-1 * ux
                                                                    + -1 * uy
                                                                    + 1 * uz)
                                                            * (-1 * ux
                                                               + -1 * uy
                                                               + 1 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 22) = (1-omega)*df_view(i.x(), i.y(), i.z(), 22)
                                                        +
                                                         0.004630f * rho *
                                                         (1.f
                                                          + 3.f * (1 * ux
                                                                   + 1 * uy
                                                                   + -1 * uz)
                                                          + 4.5f * (1 * ux
                                                                    + 1 * uy
                                                                    + -1 * uz)
                                                            * (1 * ux
                                                               + 1 * uy
                                                               + -1 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 23) = (1-omega)*df_view(i.x(), i.y(), i.z(), 23)
                                                        +
                                                         0.004630f * rho *
                                                         (1.f
                                                          + 3.f * (1 * ux
                                                                   + -1 * uy
                                                                   + -1 * uz)
                                                          + 4.5f * (1 * ux
                                                                    + -1 * uy
                                                                    + -1 * uz)
                                                            * (1 * ux
                                                               + -1 * uy
                                                               + -1 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 24) = (1-omega)*df_view(i.x(), i.y(), i.z(), 24)
                                                        +
                                                         0.004630f * rho *
                                                         (1.f
                                                          + 3.f * (-1 * ux
                                                                   + 1 * uy
                                                                   + 1 * uz)
                                                          + 4.5f * (-1 * ux
                                                                    + 1 * uy
                                                                    + 1 * uz)
                                                            * (-1 * ux
                                                               + 1 * uy
                                                               + 1 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 25) = (1-omega)*df_view(i.x(), i.y(), i.z(), 25)
                                                        +
                                                         0.004630f * rho *
                                                         (1.f
                                                          + 3.f * (-1 * ux
                                                                   + -1 * uy
                                                                   + -1 * uz)
                                                          + 4.5f * (-1 * ux
                                                                    + -1 * uy
                                                                    + -1 * uz)
                                                            * (-1 * ux
                                                               + -1 * uy
                                                               + -1 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;
                df_post_view(i.x(), i.y(), i.z(), 26) = (1-omega)*df_view(i.x(), i.y(), i.z(), 26)
                                                        +
                                                         0.004630f * rho *
                                                         (1.f
                                                          + 3.f * (1 * ux
                                                                   + 1 * uy
                                                                   + 1 * uz)
                                                          + 4.5f * (1 * ux
                                                                    + 1 * uy
                                                                    + 1 * uz)
                                                            * (1 * ux
                                                               + 1 * uy
                                                               + 1 * uz)
                                                          - 1.5f * (ux *
                                                                    ux
                                                                    + uy *
                                                                      uy
                                                                    + uz *
                                                                      uz)
                                                         ) * omega;


            }
            }
        };


        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, coll);
    }

};

#endif
