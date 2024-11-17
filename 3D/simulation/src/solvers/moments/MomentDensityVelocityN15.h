//
// Created by stloufra on 10/30/23.
//

#ifndef MOMENTDENSITYVELOCITYN15_H
#define MOMENTDENSITYVELOCITYN15_H

#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct MomentDensityVelocityN15
{
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void momentUpdate(LBMDataPointer &Data, LBMConstantsPointer &Constants) {


        auto rho_view = Data->rho.getView();

        auto ux_view = Data->ux.getView();
        auto uy_view = Data->uy.getView();
        auto uz_view = Data->uz.getView();
        
        auto df_view = Data->df.getView();
        auto mesh_view = Data->meshFluid.getView();


        auto moments = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i )mutable
        {
            if (mesh_view(i.x(), i.y(), i.z()) != 0) {

                RealType df0 = df_view(i.x(), i.y(), i.z(), 0);
                RealType df1 = df_view(i.x(), i.y(), i.z(), 1);
                RealType df2 = df_view(i.x(), i.y(), i.z(), 2);
                RealType df3 = df_view(i.x(), i.y(), i.z(), 3);
                RealType df4 = df_view(i.x(), i.y(), i.z(), 4);
                RealType df5 = df_view(i.x(), i.y(), i.z(), 5);
                RealType df6 = df_view(i.x(), i.y(), i.z(), 6);
                RealType df7 = df_view(i.x(), i.y(), i.z(), 7);
                RealType df8 = df_view(i.x(), i.y(), i.z(), 8);
                RealType df9 = df_view(i.x(), i.y(), i.z(), 9);
                RealType df10 = df_view(i.x(), i.y(), i.z(), 10);
                RealType df11 = df_view(i.x(), i.y(), i.z(), 11);
                RealType df12 = df_view(i.x(), i.y(), i.z(), 12);
                RealType df13 = df_view(i.x(), i.y(), i.z(), 13);
                RealType df14 = df_view(i.x(), i.y(), i.z(), 14);
                
                auto rho =  df0 +
                                df1 +
                                df2 +
                                df3 +
                                df4 +
                                df5 +
                                df6 +
                                df7 +
                                df8 +
                                df9 +
                                df10 +
                                df11 +
                                df12 +
                                df13 +
                                df14;

                rho_view(i.x(), i.y(), i.z()) = rho;

                ux_view(i.x(), i.y(), i.z()) = (df1
                                                  - df2
                                                  - df7
                                                  + df8
                                                  - df9
                                                  + df10
                                                  + df11
                                                  - df12
                                                  - df13
                                                  + df14) / rho;


                uy_view(i.x(), i.y(), i.z()) = (-df5
                                                  + df6
                                                  + df7
                                                  - df8
                                                  - df9
                                                  + df10
                                                  - df11
                                                  + df12
                                                  - df13
                                                  + df14) / rho;

                uz_view(i.x(), i.y(), i.z()) = (-df3
                                                  + df4
                                                  - df7
                                                  + df8
                                                  + df9
                                                  - df10
                                                  - df11
                                                  + df12
                                                  - df13
                                                  + df14) / rho;
            }
        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, moments);


    }

};

#endif
