#ifndef SOLVER_H
#define SOLVER_H

#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>
#include <omp.h>
#include <TNL/Timer.h>
#include <TNL/Logger.h>
#include <TNL/Containers/StaticArray.h>
#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>
#include "../postprocesors/outputerVTK.h"

using namespace TNL;
using namespace TNL::Algorithms;

#include "../geometry/geometryMesherBoundary.h"

#include "../traits/LBMTraits.h"

template<typename MODEL, typename MODELDATA> //TODO implement
class Solver {
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;



public:

    Solver() = delete;

    Solver(LBMConstantsPointer &Constants_,
           LBMDataPointer &Data_) {

        Constants = Constants_;
        Data = Data_;
    }

    void initializeSimulation(Vector& initial_velocity, bool verbose) {
        Constants->plot_every_it = std::ceil(Constants->plot_every / Constants->Ct);
        Constants->err_every_it = std::ceil(Constants->err_every / Constants->Ct);
        Constants->iterations = std::ceil(Constants->time / Constants->Ct);

        if (verbose) {
            std::cout << "\nPlotting every " << Constants->plot_every_it << " iterations.\n";
            std::cout << "\nCalculation will run for " << Constants->iterations << " iterations.\n";
        }

        Constants->Nvel = MODELDATA::numberOfDiscreteVelocities;

        Data->rho.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);
        Data->p.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);
        Data->u.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);
        Data->rho_out.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);
        Data->p_out.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);
        Data->u_out.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);
        Data->u_error.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);

        Data->df.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int, MODELDATA::numberOfDiscreteVelocities);
        Data->df_post.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int, MODELDATA::numberOfDiscreteVelocities);

        nonDimensionalizeBoundary(verbose);

        MODEL::initializeSim( Data, Constants, initial_velocity);

        std::cout << "Initialization of simulation done." << std::endl;
    }

    void nonDimensionalizeBoundary(bool verbose) {

        auto inlet_view = Data->meshBoundaryInlet.getView();
        auto outlet_view = Data->meshBoundaryOutlet.getView();

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

        if(false) { //TODO
            std::ofstream file("results/debug.txt");
            for (int i = 0; i < Constants->inlet_num; i++) {
                file << "x: " << inlet_view(i).velocity(0) << " ,y: " << inlet_view(i).velocity(1) << " ,z: "
                     << inlet_view(i).velocity(2) << " , i: " << i << std::endl;
            }
            file.close();
        }

        if (verbose) {
            std::cout << "Inlet non-dimensionalized.\n";
        }

        auto nonDimOutlet = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {
            outlet_view(i.x()).density =
                    outlet_view(i.x()).density / Cm * Cl * Cl * Cl;
        };

        parallelFor<DeviceType>(0, Constants->outlet_num, nonDimOutlet);

        if (verbose) {
            std::cout << "Outlet non-dimensionalized.\n";
        }

    }

    void runSimulation() {

        timer_loop.start();


        std::cout << "Got here.\n";



        int k = 0;
        while(k<Constants -> iterations) //err>=10e-4)
        {


            timer_collision.start();
            MODEL::collision(Data, Constants);
            timer_collision.stop();

            timer_streaming.start();
            MODEL::streaming(Data, Constants);
            timer_streaming.stop();

            timer_bounceback.start();
            MODEL::bounceBack(Data, Constants);
            timer_bounceback.stop();

            timer_momentsUpdate.start();
            MODEL::momentUpdate(Data, Constants);
            timer_momentsUpdate.stop();


            if(k%Constants -> err_every_it==0 && k!=0)
            {
                timer_err.start();
                MODEL::errorEvaluation(Data, Constants);
                //TODO ux_center=%e uy_center=%e rho_center=%e
                printf("\n err=%e k=%d\n",Constants->err, k);
                if (std::isnan(Constants->err))
                {
                    std::cout << "\n Error is NaN, breaking out.\n";
                    break;
                }
                timer_err.stop();
            }

            if(k%Constants -> plot_every_it==0)
            {
                timer_momentsUpdate.start();
                MODEL::pressureUpdate(Data, Constants);
                timer_momentsUpdate.stop();

                timer_output.start();
                outputerVTK::variablesLatticeVTK(Data, Constants, k/Constants -> plot_every_it, 1);
                timer_output.stop();
            }

            k++;

        }


    };

    void convertToLattice(bool verbose) {
        // prenasobim lattice jednotky a dostanu fyzikalni
        // vydelim fyz a dostanu lattice
        if (Constants->U_lb > 0.1) {
            std::cout << "\n Lattice speed should not be higher than 0.1 due to stability close to speed of sound.\n";
        }

        assert(Constants->U_lb < Constants->cs);

        Constants->L_fyz = abs(Constants->BBmaxx - Constants->BBminx) * Constants->conversion_factor_fyz;

        Constants->Cl = Constants->L_fyz / Constants->dimX_int;

        RealType U_fyz;
        RealType U_mag;

        /*TODO std::ofstream file("results/debug.txt");
        for (int i = 0; i < Constants->inlet_num; i++) {
            boundaryConditionInlet BC = Data->meshBoundaryInlet[i];
            U_mag = TNL::l2Norm(
                    BC.velocity);// sqrt(BC.velocity.x * BC.velocity.x + BC.velocity.y * BC.velocity.y +BC.velocity.z * BC.velocity.z);

                      U_fyz = std::max(U_mag, U_fyz);

            file <<"x: " <<BC.velocity.x()<<" ,y: " << BC.velocity.y()<< " ,z: "   <<BC.velocity.z() << " ,normUfyz: " << U_fyz <<  std::endl;
        }

        U_fyz = std::max(U_fyz,Constants -> u_guess_fyz);
        file.close();
         */
        U_fyz = Constants -> u_guess_fyz;
       // if (verbose) {
        //    std::cout << "\n- Maximal velocity is " << U_fyz << ".\n";
        //}

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


    LBMConstantsPointer Constants;
    LBMDataPointer Data;


    Timer timer_loop;
    Timer timer_collision;
    Timer timer_streaming;
    Timer timer_bounceback;
    Timer timer_momentsUpdate;
    Timer timer_err;
    Timer timer_output;

};

#endif //SOLVER_H
