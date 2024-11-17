//
// Created by stloufra on 10/30/23.
//

#ifndef STREAMINGABPULL_H
#define STREAMINGABPULL_H

#include "../../traits/LBMTraits.h"

template <typename MODELDATA>
struct StreamingABpull
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

            if (mesh_val > 0 || mesh_val == -2 || mesh_val == -3) // pull schema
            {
                for (int vel = 0; vel < Nvel; vel++)
                {
                    int id;
                    int jd;
                    int kd;
                    id = i.x() - MD.c[vel][0];
                    jd = i.y() - MD.c[vel][1];
                    kd = i.z() - MD.c[vel][2];


                    auto mesh_neigh_val = mesh_view(id, jd, kd);

                    if (mesh_neigh_val > 0 || mesh_view(id, jd,  kd) == -2 || mesh_neigh_val == -3)
                    {
                        df_view(i.x(), i.y(), i.z(), vel) = df_post_view(id, jd, kd, vel);
                    }
                    else if (mesh_neigh_val == -1)
                    {
                        df_view(i.x(), i.y(), i.z(), vel) = df_view(id, jd, kd, vel);
                    }
                }

                if(mesh_val == -3){ //push
                  for (int vel = 0; vel < Nvel; vel++){
                    int id;
                    int jd;
                    int kd;
                    id = i.x() + MD.c[vel][0];
                    jd = i.y() + MD.c[vel][1];
                    kd = i.z() + MD.c[vel][2];

                    df_view(id, jd, kd, vel) = df_post_view(i.x(), i.y(), i.z(), vel);

                    }
                }
            }
        };

        TNL::Containers::StaticArray < 3, int > begin{0, 0, 0};
        TNL::Containers::StaticArray < 3, int > end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};

        //DeviceType::LaunchConfiguration A;
        //A.blockSize = (4, 4, 4);
        parallelFor<DeviceType>(begin, end, stream);
    }
};

#endif
