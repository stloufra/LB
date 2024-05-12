//
// Created by stloufra on 10/30/23.
//

#ifndef OMEGANO_H
#define OMEGANO_H

#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct OmegaNo {


    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using TensorType = LBMTraits::TensorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void omega(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        auto omega_view = Data->omega.getView();
        auto mesh_view = Data->meshFluid.getView();
        const auto tau = Constants->tau;

        MODELDATA MD;


        auto omegaFunc = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i ) mutable
        {
            if (mesh_view(i.x(), i.y(), i.z()) != 0) {

                omega_view(i.x(), i.y(), i.z()) = 1.f / tau;
            }
        };


        TNL::Containers::StaticArray<3, int> begin1{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end1{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin1, end1, omegaFunc);
    }


};

#endif
