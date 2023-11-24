//
// Created by stloufra on 10/30/23.
//

#ifndef STREAMINGAB_H
#define STREAMINGAB_H

#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct StreamingAB
{
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void streaming(LBMDataPointer &Data, LBMConstantsPointer &Constants) {


        using RealType = LBMTraits::RealType;
        using DeviceType = LBMTraits::DeviceType;
        using VectorType = LBMTraits::VectorType;
        using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
        using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

        auto mesh_ConstView = Data->meshFluid.getConstView();
        auto df_post_ConstView = Data->df_post.getConstView();

        auto df_view = Data->df.getView();

        const auto Nvel = Constants->Nvel;
        const auto dimX_int = Constants->dimX_int;
        const auto dimY_int = Constants->dimY_int;
        const auto dimZ_int = Constants->dimZ_int;

        MODELDATA MD;

        auto stream = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i  ) mutable
        {
            if (mesh_ConstView(i.x(), i.y(), i.z()) != 0) { //TODO

                int id;
                int jd;
                int kd;


                for (int vel = 0; vel < Nvel; vel++) {

                    id = i.x() + MD.c[vel][0];
                    jd = i.y() + MD.c[vel][1];
                    kd = i.z() + MD.c[vel][2];


                    //if (id >= 0 && jd >= 0 && kd >= 0 && id < dimX_int && jd < dimY_int &&
                    //kd < dimZ_int) {
                        if (mesh_ConstView(id, jd, kd) != 0) {
                            df_view(id, jd, kd, vel) = df_post_ConstView(i.x(), i.y(), i.z(), vel);
                        }
                    //}
                }
            }

        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{dimX_int, dimY_int, dimZ_int};
        parallelFor<DeviceType>(begin, end, stream);
    }

};

#endif
