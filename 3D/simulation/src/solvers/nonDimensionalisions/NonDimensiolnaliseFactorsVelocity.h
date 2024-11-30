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

        assert(Constants->U_lb < Constants->cs && "Speed should not be greater than speed of sound!");

        Constants->L_fyz = abs(Constants->BBmaxx - Constants->BBminx) * Constants->conversion_factor_fyz;
        RealType L_fyz_y = abs(Constants->BBmaxy - Constants->BBminy) * Constants->conversion_factor_fyz;
        RealType L_fyz_z = abs(Constants->BBmaxz - Constants->BBminz) * Constants->conversion_factor_fyz;

        std::cout << "L in X direction is: " <<  Constants->L_fyz << std::endl;
        std::cout << "L in Y direction is: " << L_fyz_y << std::endl;
        std::cout << "L in Z direction is: " << L_fyz_z << std::endl;


        Constants->Cl = Constants->L_fyz / Constants->dimX_int;

        RealType U_fyz;

        U_fyz = Constants -> u_guess_fyz; //TODO: think of better way to do this

        Constants->Cu = U_fyz / Constants->U_lb;
        Constants->Cu_inverse = 1 / Constants->Cu;
        Constants->Ct = Constants->Cl / Constants->Cu;
        Constants->Crho = Constants->rho_fyz;
        Constants->Cm = Constants->Crho * Constants->Cl * Constants->Cl * Constants->Cl;
        Constants->Cpressure = Constants->Cm/Constants->Cl/Constants->Ct/Constants->Ct;

        Constants->Re = Constants->U_inf * Constants->L_fyz / Constants->ny_fyz;
        RealType Re_y = Constants->U_inf * L_fyz_y / Constants->ny_fyz;
        RealType Re_z = Constants->U_inf * L_fyz_z / Constants->ny_fyz;

        Constants->U_inf = Constants->U_inf * Constants->Cu_inverse;

        std::cout << "\n- U_inf  is " << Constants->U_inf << "\n";

        std::cout << "\n- Re (for L in X direction) is " << Constants->Re << "\n";
        std::cout << "\n- Re (for L in Y direction) is " << Re_y << "\n";
        std::cout << "\n- Re (for L in Z direction) is " << Re_z << "\n";

        Constants->ny = Constants->ny_fyz * Constants->Ct / Constants->Cl / Constants->Cl;

        std::cout << "\n- Nu  is " << Constants->ny << "\n";

        Constants->tau = Constants->ny / 3.f + 0.5f;

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

        Constants->omega = 1.f / Constants->tau;

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

    static void nonDimensionalizePeriodicDP(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        auto periodic_view = Data->meshBoundaryPeriodicDP.getView();

        auto Cu_inverse = Constants->Cu_inverse;
        auto Cm= Constants->Cm;
        auto Cl = Constants->Cl;

        auto cs2 = Constants ->cs2;

        auto Cpressure = Constants->Cpressure;

        //printf("Nondimesionalized delta Rho - %f", periodic_view(0).DeltaRho);

        auto nonDimPerDP = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &nod  ) mutable
        {
            periodic_view(nod.x()).DeltaRho =
                    periodic_view(nod.x()).DeltaRho / cs2  / Cpressure;
        };

        parallelFor<DeviceType>(0, Constants->periodicDP_num, nonDimPerDP);

        //printf("Nondimesionalized delta Rho - %f", periodic_view(0).DeltaRho);

    }

};

#endif
