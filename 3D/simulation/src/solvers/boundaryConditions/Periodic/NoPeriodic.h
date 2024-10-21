//
// Created by stloufra on 10/30/23.
//

#ifndef NOPERIODIC_H
#define NOPERIODIC_H

#include "../../../traits/LBMTraits.h"

template<typename MODELDATA>
struct NoPeriodic {
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void periodic(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

    }
};

#endif
