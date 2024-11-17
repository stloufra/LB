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

        auto uxStatArr_view = Data->ux.getStorageArray().getConstView();
        auto uyStatArr_view = Data->uy.getStorageArray().getConstView();
        auto uzStatArr_view = Data->uz.getStorageArray().getConstView();
        
        auto uxErrorStatArr_view = Data->ux_error.getStorageArray().getConstView();
        auto uyErrorStatArr_view = Data->uy_error.getStorageArray().getConstView();
        auto uzErrorStatArr_view = Data->uz_error.getStorageArray().getConstView();


        auto mesh_view = Data->meshFluid.getStorageArray().getConstView();

        auto ux_view = Data->ux.getView();
        auto uy_view = Data->uy.getView();
        auto uz_view = Data->uz.getView();


        auto uxError_view = Data->ux_error.getView();
        auto uyError_view = Data->uy_error.getView();
        auto uzError_view = Data->uz_error.getView();


        RealType err1, err2;
        err1 = 0.0f;
        err2 = 0.0f;

        auto fetch1 = [=]
        __cuda_callable__(int
        i ) ->RealType
        {
            //dui^2

            auto ux = uxStatArr_view[i];
            auto uy = uyStatArr_view[i];
            auto uz = uzStatArr_view[i];

            auto uxE = uxErrorStatArr_view[i];
            auto uyE = uyErrorStatArr_view[i];
            auto uzE = uzErrorStatArr_view[i];

            if (mesh_view[i] == 1) {
                const RealType &add1 = (ux - uxE) *
                                       (ux - uxE) +
                                       (uy - uyE) *
                                       (uy - uyE) +
                                       (uz - uzE) *
                                       (uz - uzE);

                return sqrt(add1);
            }
            return 0;
        };

        auto fetch2 = [=]
        __cuda_callable__(int
        i ) ->RealType
        {
            //ui^2

            auto ux = uxStatArr_view[i];
            auto uy = uyStatArr_view[i];
            auto uz = uzStatArr_view[i];

            if (mesh_view[i] == 1) {
                const RealType &add2 =
                        ux * ux +
                        uy * uy +
                        uz * uz;

                return sqrt(add2);
            }
            return 0;
        };

        err1 = sqrt(reduce<DeviceType, int>(0,
                                            uxStatArr_view.getSize(),
                                            fetch1,
                                            TNL::Plus(),
                                            0.0));

        err2 = sqrt(reduce<DeviceType, int>(0,
                                            uxStatArr_view.getSize(),
                                            fetch2,
                                            TNL::Plus(),
                                            0.0));
        Constants->err = err1 / err2;

        auto copy = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i  ) mutable
        {
            uxError_view(i.x(), i.y(), i.z()) = ux_view(i.x(), i.y(), i.z());
            uyError_view(i.x(), i.y(), i.z()) = uy_view(i.x(), i.y(), i.z());
            uzError_view(i.x(), i.y(), i.z()) = uz_view(i.x(), i.y(), i.z());

        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, copy);
    }

};

#endif
