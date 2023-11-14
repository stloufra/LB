//
// Created by stloufra on 10/30/23.
//

#ifndef ERRORQUADRATIC_H
#define ERRORQUADRATIC_H

#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct ErrorQuadratic
{
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void errorEvaluation(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        auto uStatArr_view = Data->u.getStorageArray().getConstView();
        auto uErrorStatArr_view = Data->u_error.getStorageArray().getConstView();
        auto mesh_view = Data->meshFluid.getStorageArray().getConstView();

        auto u_view = Data->u.getView();
        auto uError_view = Data->u_error.getView();


        RealType err1, err2;
        err1 = 0.0f;
        err2 = 0.0f;

        auto fetch1 = [=]
        __cuda_callable__(int
        i ) ->RealType
        {
            //dui^2
            if (mesh_view[i] == 1) {
                const RealType &add1 = (uStatArr_view[i].x() - uErrorStatArr_view[i].x()) *
                                       (uStatArr_view[i].x() - uErrorStatArr_view[i].x()) +
                                       (uStatArr_view[i].y() - uErrorStatArr_view[i].y()) *
                                       (uStatArr_view[i].y() - uErrorStatArr_view[i].y()) +
                                       (uStatArr_view[i].y() - uErrorStatArr_view[i].z()) *
                                       (uStatArr_view[i].y() - uErrorStatArr_view[i].z());

                return sqrt(add1);
            }
            return 0;
        };

        auto fetch2 = [=]
        __cuda_callable__(int
        i ) ->RealType
        {
            //ui^2
            if (mesh_view[i] == 1) {
                const RealType &add2 =
                        uStatArr_view[i].x() * uStatArr_view[i].x() +
                        uStatArr_view[i].y() * uStatArr_view[i].y() +
                        uStatArr_view[i].z() * uStatArr_view[i].z();

                return sqrt(add2);
            }
            return 0;
        };

        err1 = sqrt(reduce<DeviceType, int>(0,
                                            uStatArr_view.getSize(),
                                            fetch1,
                                            TNL::Plus(),
                                            0.0));

        err2 = sqrt(reduce<DeviceType, int>(0,
                                            uStatArr_view.getSize(),
                                            fetch2,
                                            TNL::Plus(),
                                            0.0));
        Constants->err = err1 / err2;

        auto copy = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i  ) mutable
        {
            uError_view(i.x(), i.y(), i.z()) = u_view(i.x(), i.y(), i.z());

        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, copy);
    }

};

#endif
