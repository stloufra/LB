//
// Created by stloufra on 10/30/23.
//

#ifndef NOSYMMETRY_H
#define NOSYMMETRY_H

#include "../../../traits/LBMTraits.h"

template<typename MODELDATA>
struct NoSymmetry {
    using DeviceType = LBMTraits::DeviceType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void symmetry(LBMDataPointer &Data, LBMConstantsPointer &Constants){

    };

};

#endif
