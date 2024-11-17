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

        auto ux_view = Data->ux.getView();
        auto uy_view = Data->uy.getView();
        auto uz_view = Data->uz.getView();

        auto rhoTimeAvg_view = Data->rhoTimeAvg.getView();

        auto uxTimeAvg_view = Data->uxTimeAvg.getView();
        auto uyTimeAvg_view = Data->uyTimeAvg.getView();
        auto uzTimeAvg_view = Data->uzTimeAvg.getView();


        auto timeAverage = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i)
        mutable
        {

            rhoTimeAvg_view(i.x(), i.y(), i.z()) += rho_view(i.x(), i.y(), i.z());

            uxTimeAvg_view(i.x(), i.y(), i.z()) += ux_view(i.x(), i.y(), i.z());
            uyTimeAvg_view(i.x(), i.y(), i.z()) += uy_view(i.x(), i.y(), i.z());
            uzTimeAvg_view(i.x(), i.y(), i.z()) += uz_view(i.x(), i.y(), i.z());

        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, timeAverage);

         Constants-> timeAveraged = false;

        Constants-> TimeAvgCounter += 1;

    }

    static void momentAvg(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        auto rhoTimeAvg_view = Data->rhoTimeAvg.getView();

        auto uxTimeAvg_view = Data->uxTimeAvg.getView();
        auto uyTimeAvg_view = Data->uyTimeAvg.getView();
        auto uzTimeAvg_view = Data->uzTimeAvg.getView();

        RealType TAC = static_cast<RealType>(Constants->iterationsMomentAvg);

        printf("Counter %d" , Constants->TimeAvgCounter);

        auto timeAverage = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i)
        mutable
        {

            rhoTimeAvg_view(i.x(), i.y(), i.z()) /= TAC;

            uxTimeAvg_view(i.x(), i.y(), i.z())/= TAC;
            uyTimeAvg_view(i.x(), i.y(), i.z())/= TAC;
            uzTimeAvg_view(i.x(), i.y(), i.z())/= TAC;

        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, timeAverage);

         Constants-> timeAveraged = true;
    }

    static void momentAvgDelete(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        auto rhoTimeAvg_view = Data->rhoTimeAvg.getView();

        auto uxTimeAvg_view = Data->uxTimeAvg.getView();
        auto uyTimeAvg_view = Data->uyTimeAvg.getView();
        auto uzTimeAvg_view = Data->uzTimeAvg.getView();

        auto timeAveraged = Constants->timeAveraged;

        auto timeAverage = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i)
        mutable
        {

            rhoTimeAvg_view(i.x(), i.y(), i.z()) = 0;

            uxTimeAvg_view(i.x(), i.y(), i.z()) = 0;
            uyTimeAvg_view(i.x(), i.y(), i.z()) = 0;
            uzTimeAvg_view(i.x(), i.y(), i.z()) = 0;

        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, timeAverage);

        Constants-> TimeAvgCounter == 0; //for another one
        }
};

#endif
