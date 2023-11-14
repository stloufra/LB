//
// Created by stloufra on 10/30/23.
//

#ifndef MOMENTDENSITYVELOCITYN27_H
#define MOMENTDENSITYVELOCITYN27_H

#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct MomentDensityVelocityN27
{
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void momentUpdate(LBMDataPointer &Data, LBMConstantsPointer &Constants) {


        auto rho_view = Data->rho.getView();
        auto u_view = Data->u.getView();
        auto df_view = Data->df.getView();
        auto mesh_view = Data->meshFluid.getView();


        auto moments = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i )mutable
        {
            if (mesh_view(i.x(), i.y(), i.z()) != 0) {
                rho_view(i.x(), i.y(), i.z()) = df_view(i.x(), i.y(), i.z(), 0) +
                                                df_view(i.x(), i.y(), i.z(), 1) +
                                                df_view(i.x(), i.y(), i.z(), 2) +
                                                df_view(i.x(), i.y(), i.z(), 3) +
                                                df_view(i.x(), i.y(), i.z(), 4) +
                                                df_view(i.x(), i.y(), i.z(), 5) +
                                                df_view(i.x(), i.y(), i.z(), 6) +
                                                df_view(i.x(), i.y(), i.z(), 7) +
                                                df_view(i.x(), i.y(), i.z(), 8) +
                                                df_view(i.x(), i.y(), i.z(), 9) +
                                                df_view(i.x(), i.y(), i.z(), 10) +
                                                df_view(i.x(), i.y(), i.z(), 11) +
                                                df_view(i.x(), i.y(), i.z(), 12) +
                                                df_view(i.x(), i.y(), i.z(), 13) +
                                                df_view(i.x(), i.y(), i.z(), 14) +
                                                df_view(i.x(), i.y(), i.z(), 15) +
                                                df_view(i.x(), i.y(), i.z(), 16) +
                                                df_view(i.x(), i.y(), i.z(), 17) +
                                                df_view(i.x(), i.y(), i.z(), 18) +
                                                df_view(i.x(), i.y(), i.z(), 19) +
                                                df_view(i.x(), i.y(), i.z(), 20) +
                                                df_view(i.x(), i.y(), i.z(), 21) +
                                                df_view(i.x(), i.y(), i.z(), 22) +
                                                df_view(i.x(), i.y(), i.z(), 23) +
                                                df_view(i.x(), i.y(), i.z(), 24) +
                                                df_view(i.x(), i.y(), i.z(), 25) +
                                                df_view(i.x(), i.y(), i.z(), 26);

                u_view(i.x(), i.y(), i.z())(0) = (df_view(i.x(), i.y(), i.z(), 1)
                                                  - df_view(i.x(), i.y(), i.z(), 2)
                                                  + df_view(i.x(), i.y(), i.z(), 7)
                                                  - df_view(i.x(), i.y(), i.z(), 8)
                                                  + df_view(i.x(), i.y(), i.z(), 9)
                                                  - df_view(i.x(), i.y(), i.z(), 10)
                                                  - df_view(i.x(), i.y(), i.z(), 11)
                                                  + df_view(i.x(), i.y(), i.z(), 12)
                                                  - df_view(i.x(), i.y(), i.z(), 15)
                                                  + df_view(i.x(), i.y(), i.z(), 16)
                                                  - df_view(i.x(), i.y(), i.z(), 19)
                                                  + df_view(i.x(), i.y(), i.z(), 20)
                                                  - df_view(i.x(), i.y(), i.z(), 21)
                                                  + df_view(i.x(), i.y(), i.z(), 22)
                                                  + df_view(i.x(), i.y(), i.z(), 23)
                                                  - df_view(i.x(), i.y(), i.z(), 24)
                                                  - df_view(i.x(), i.y(), i.z(), 25)
                                                  + df_view(i.x(), i.y(), i.z(), 26)) / rho_view(i.x(), i.y(), i.z());


                u_view(i.x(), i.y(), i.z())(1) = (-df_view(i.x(), i.y(), i.z(), 5)
                                                  + df_view(i.x(), i.y(), i.z(), 6)
                                                  - df_view(i.x(), i.y(), i.z(), 11)
                                                  + df_view(i.x(), i.y(), i.z(), 12)
                                                  + df_view(i.x(), i.y(), i.z(), 13)
                                                  - df_view(i.x(), i.y(), i.z(), 14)
                                                  + df_view(i.x(), i.y(), i.z(), 15)
                                                  - df_view(i.x(), i.y(), i.z(), 16)
                                                  + df_view(i.x(), i.y(), i.z(), 17)
                                                  - df_view(i.x(), i.y(), i.z(), 18)
                                                  + df_view(i.x(), i.y(), i.z(), 19)
                                                  - df_view(i.x(), i.y(), i.z(), 20)
                                                  - df_view(i.x(), i.y(), i.z(), 21)
                                                  + df_view(i.x(), i.y(), i.z(), 22)
                                                  - df_view(i.x(), i.y(), i.z(), 23)
                                                  + df_view(i.x(), i.y(), i.z(), 24)
                                                  - df_view(i.x(), i.y(), i.z(), 25)
                                                  + df_view(i.x(), i.y(), i.z(), 26)) / rho_view(i.x(), i.y(), i.z());

                u_view(i.x(), i.y(), i.z())(2) = (-df_view(i.x(), i.y(), i.z(), 3)
                                                  + df_view(i.x(), i.y(), i.z(), 4)
                                                  - df_view(i.x(), i.y(), i.z(), 7)
                                                  + df_view(i.x(), i.y(), i.z(), 8)
                                                  + df_view(i.x(), i.y(), i.z(), 9)
                                                  - df_view(i.x(), i.y(), i.z(), 10)
                                                  - df_view(i.x(), i.y(), i.z(), 13)
                                                  + df_view(i.x(), i.y(), i.z(), 14)
                                                  + df_view(i.x(), i.y(), i.z(), 17)
                                                  - df_view(i.x(), i.y(), i.z(), 18)
                                                  - df_view(i.x(), i.y(), i.z(), 19)
                                                  + df_view(i.x(), i.y(), i.z(), 20)
                                                  + df_view(i.x(), i.y(), i.z(), 21)
                                                  - df_view(i.x(), i.y(), i.z(), 22)
                                                  - df_view(i.x(), i.y(), i.z(), 23)
                                                  + df_view(i.x(), i.y(), i.z(), 24)
                                                  - df_view(i.x(), i.y(), i.z(), 25)
                                                  + df_view(i.x(), i.y(), i.z(), 26)) / rho_view(i.x(), i.y(), i.z());
            }
        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, moments);


    }

};

#endif
