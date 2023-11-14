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

        auto mesh_view = Data->meshFluid.getView();
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();

        const auto Nvel = Constants->Nvel;
        auto dimX_int = Constants->dimX_int;
        auto dimY_int = Constants->dimY_int;
        auto dimZ_int = Constants->dimZ_int;

        MODELDATA MD;

        auto stream = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i  ) mutable
        {
            if (mesh_view(i.x(), i.y(), i.z()) != 0) { //TODO

                for (int vel = 0; vel < Nvel; vel++) {
                    int id;
                    int jd;
                    int kd;
                    id = i.x() + MD.c[vel][0];
                    jd = i.y() + MD.c[vel][1];
                    kd = i.z() + MD.c[vel][2];


                    if (id >= 0 && jd >= 0 && kd >= 0 && id < dimX_int && jd < dimY_int &&
                        kd < dimZ_int) {
                        if (mesh_view(id, jd, kd) != 0) {
                            df_view(id, jd, kd, vel) = df_post_view(i.x(), i.y(), i.z(), vel);
                        }
                    }
                }
            }

        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, stream);
    }

};

#endif
