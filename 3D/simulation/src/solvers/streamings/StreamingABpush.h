//
// Created by stloufra on 10/30/23.
//

#ifndef STREAMINGABPUSH_H
#define STREAMINGABPUSH_H

#include "../../traits/LBMTraits.h"

template <typename MODELDATA>
struct StreamingABpush
{
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void streaming(LBMDataPointer& Data, LBMConstantsPointer& Constants)
    {
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
                const TNL::Containers::StaticArray<3, int>& i) mutable{

            auto mesh_val = mesh_view(i.x(), i.y(), i.z());

            if (mesh_val != 0 && mesh_val != -1) // push schema
            {
                for (int vel = 0; vel < Nvel; vel++)
                {
                    int id;
                    int jd;
                    int kd;
                    id = i.x() + MD.c[vel][0];
                    jd = i.y() + MD.c[vel][1];
                    kd = i.z() + MD.c[vel][2];

                    df_view(id, jd, kd, vel) = df_post_view(i.x(), i.y(), i.z(), vel);

                }
            }
            else if (mesh_val == -1)
            {
                for (int vel = 0; vel < Nvel; vel++)
                {
                    int id;
                    int jd;
                    int kd;
                    id = i.x() + MD.c[vel][0];
                    jd = i.y() + MD.c[vel][1];
                    kd = i.z() + MD.c[vel][2];

                df_view(id, jd, kd, vel) = df_view(i.x(), i.y(), i.z(), vel);
                }
            }
        };

        TNL::Containers::StaticArray < 3, int > begin{1, 1, 1};
        TNL::Containers::StaticArray < 3, int > end{Constants->dimX_int-1, Constants->dimY_int-1, Constants->dimZ_int-1};
        //DeviceType::LaunchConfiguration A;
        //A.blockSize = (128, 128, 128);
        parallelFor<DeviceType>(begin, end, stream);
        //std::cout <<A.blockSize.x <<", "<<A.blockSize.y <<", "<<A.blockSize.z <<", " << std::endl ;
    }
};

#endif
