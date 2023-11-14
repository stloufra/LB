//
// Created by stloufra on 10/30/23.
//

#ifndef MOMENTPRESSURE_H
#define MOMENTPRESSURE_H

#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct MomentPressure
{
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void momentUpdate(LBMDataPointer &Data, LBMConstantsPointer &Constants) {
        auto p_view = Data->p.getView();
        auto rho_view = Data->rho.getView();
        auto mesh_view = Data->meshFluid.getView();
        auto cs2 = Constants->cs2;


        auto pressure = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i)
        mutable
        {

            p_view(i.x(), i.y(), i.z()) = rho_view(i.x(), i.y(), i.z()) * cs2;

        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, pressure);
    }
};

#endif
