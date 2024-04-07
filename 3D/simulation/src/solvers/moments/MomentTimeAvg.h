//
// Created by stloufra on 10/30/23.
//

#ifndef MOMENTTIMEAVG_H
#define MOMENTTIMEAVG_H

#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct MomentTimeAvg
{
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void momentAdd(LBMDataPointer &Data, LBMConstantsPointer &Constants) {
        auto rho_view = Data->rho.getView();
        auto u_view = Data->u.getView();

        auto rhoTimeAvg_view = Data->rhoTimeAvg.getView();
        auto uTimeAvg_view = Data->uTimeAvg.getView();


        auto timeAverage = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i)
        mutable
        {

            rhoTimeAvg_view(i.x(), i.y(), i.z()) += rho_view(i.x(), i.y(), i.z());

            uTimeAvg_view(i.x(), i.y(), i.z()).x() += u_view(i.x(), i.y(), i.z()).x();
            uTimeAvg_view(i.x(), i.y(), i.z()).y() += u_view(i.x(), i.y(), i.z()).y();
            uTimeAvg_view(i.x(), i.y(), i.z()).z() += u_view(i.x(), i.y(), i.z()).z();

        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, timeAverage);

        Constants->TimeAvgCounter += 1;

    }

    static void momentAvg(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        auto rhoTimeAvg_view = Data->rhoTimeAvg.getView();
        auto uTimeAvg_view = Data->uTimeAvg.getView();

        auto timeAveraged = Constants->timeAveraged;

        auto timeAverage = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i)
        mutable
        {

            rhoTimeAvg_view(i.x(), i.y(), i.z()) /= Constants->TimeAvgCounter;

            uTimeAvg_view(i.x(), i.y(), i.z()).x() /= Constants->TimeAvgCounter;
            uTimeAvg_view(i.x(), i.y(), i.z()).y() /= Constants->TimeAvgCounter;
            uTimeAvg_view(i.x(), i.y(), i.z()).z() /= Constants->TimeAvgCounter;

        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, timeAverage);


        timeAveraged = true;
    }
};

#endif
