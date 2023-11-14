//
// Created by stloufra on 10/30/23.
//

#ifndef NONDIMANSIONALIZEFACTROSVELOCITY_H
#define NONDIMANSIONALIZEFACTROSVELOCITY_H

#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct NonDimensiolnaliseFactorsVelocity
{

    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void nonDimensionalize(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        // prenasobim lattice jednotky a dostanu fyzikalni
        // vydelim fyz a dostanu lattice

        if (Constants->U_lb > 0.1) {
            std::cout << "\n Lattice speed should not be higher than 0.1 due to stability close to speed of sound.\n";
        }

        assert(Constants->U_lb < Constants->cs);

        Constants->L_fyz = abs(Constants->BBmaxx - Constants->BBminx) * Constants->conversion_factor_fyz;

        Constants->Cl = Constants->L_fyz / Constants->dimX_int;

        RealType U_fyz;

        U_fyz = Constants -> u_guess_fyz;

        Constants->Cu = U_fyz / Constants->U_lb;
        Constants->Cu_inverse = 1 / Constants->Cu;
        Constants->Ct = Constants->Cl / Constants->Cu;
        Constants->Crho = Constants->rho_fyz;
        Constants->Cm = Constants->Crho * Constants->Cl * Constants->Cl * Constants->Cl;
        Constants->Cpressure = Constants->Cm/Constants->Cl/Constants->Ct/Constants->Ct;
        Constants->Re = U_fyz * Constants->L_fyz / Constants->ny_fyz;

        std::cout << "\n- Re is " << Constants->Re << "\n";

        Constants->ny = Constants->ny_fyz * Constants->Ct / Constants->Cl / Constants->Cl;

        Constants->tau = Constants->ny / 3 + 0.5f;

        if (Constants->tau < 0.51) {
            std::cout << "Tau is too small. Consider higher resolution or higher lattice speed. \n";
            std::cout << "Tau is " << Constants->tau << "\n";
        } else if (Constants->tau > 0.99) {
            std::cout << "Tau is too high. Consider lower resolution or lower lattice speed. \n";
            std::cout << "Tau is " << Constants->tau << "\n";
        } else {
            std::cout << "- Tau is " << Constants->tau << " <-(0.5;1)\n";

        }

        assert(Constants->tau > 0.5 && Constants->tau < 0.99);

        Constants->omega = 1 / Constants->tau;

        std::cout << "\nConversion undergone successfully." << "\n";

    }

    static void nonDimensionalizeInlet(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        auto inlet_view = Data->meshBoundaryInlet.getView();

        auto Cu_inverse = Constants->Cu_inverse;
        auto Cm= Constants->Cm;
        auto Cl = Constants->Cl;



        auto nonDimInlet = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {
            inlet_view(i.x()).velocity(0) = inlet_view(i.x()).velocity.x() * Cu_inverse;
            inlet_view(i.x()).velocity(1) = inlet_view(i.x()).velocity.y() * Cu_inverse;
            inlet_view(i.x()).velocity(2)=  inlet_view(i.x()).velocity.z() * Cu_inverse;

        };

        parallelFor<DeviceType>(0, Constants->inlet_num, nonDimInlet);


    }

    static void nonDimensionalizeOutlet(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        auto outlet_view = Data->meshBoundaryOutlet.getView();

        auto Cu_inverse = Constants->Cu_inverse;
        auto Cm= Constants->Cm;
        auto Cl = Constants->Cl;


        auto nonDimOutlet = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {
            outlet_view(i.x()).density =
                    outlet_view(i.x()).density / Cm * Cl * Cl * Cl;
        };

        parallelFor<DeviceType>(0, Constants->outlet_num, nonDimOutlet);

    }

};

#endif
