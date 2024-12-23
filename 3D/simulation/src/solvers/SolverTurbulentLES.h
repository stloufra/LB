#ifndef SOLVERTURBULENTLES_H
#define SOLVERTURBULENTLES_H

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

#include "./collisions/CollisionSRTTurbulent.h"
#include "./collisions/CollisionCumD3Q27Turbulent2015.h"
#include "./collisions/CollisionCumD3Q27Turbulent2017.h"
#include "./collisions/CollisionCumD3Q27TurbulentCombined.h"

#include "./initializations/InitializationEquilibriumConstVector.h"
#include "./initializations/InitializationEquilibriumVariables.h"

#include "./streamings//StreamingAB.h"

#include "boundaryConditions/Wall/BounceBackWallHalf.h"

#include "boundaryConditions/Periodic/Periodic.h"
#include "boundaryConditions/Periodic/NoPeriodic.h"

#include "boundaryConditions/Symmetry/BounceSymmetryHalf.h"
#include "boundaryConditions/Symmetry/NoSymmetry.h"

#include "boundaryConditions/Inlet/InletVelocityMovingWall.h"
#include "boundaryConditions/Inlet/InletVelocityZouHe.h"
#include "boundaryConditions/Inlet/InletVelocityEquilibrium.h"

#include "./boundaryConditions/Outlet/OutletDensityInterpolated.h"
#include "boundaryConditions/Outlet/OutletDensityEquilibrium.h"
#include "boundaryConditions/Outlet/OutletNeighbourEquilibrium.h"
#include "boundaryConditions/Outlet/OutletNeighbourEquilibriumOmega.h"
#include "boundaryConditions/Outlet/OutletNeighbourEquilibriumOmegaRF.h"

#include "./moments/MomentDensityVelocityN27.h"
#include "./moments/MomentDensityVelocityN19.h"
#include "./moments/MomentDensityVelocityN15.h"
#include "./moments/MomentPressure.h"
#include "./moments/MomentTimeAvg.h"

#include "./turbulence/OmegaLES.h"
#include "./turbulence/OmegaNo.h"

#include "./errorEvaluations/ErrorQuadratic.h"

#include "./nonDimensionalisions/NonDimensiolnaliseFactorsVelocity.h"

#include "../geometry/geometryMesherBoundary.h"
#include "../postprocesors/outputerVTK.h"
#include "../traits/LBMTraits.h"

using namespace TNL;
using namespace TNL::Algorithms;

template<   typename MODELTYPE,
            typename INITIALIZATIONTYPE,
            typename COLLISIONTYPE,
            typename STREAMINGTYPE,
            typename BOUNCEBACKWALLTYPE,
            typename SYMMETRYTYPE,
            typename PERIODICTYPE,
            typename INLETTYPE,
            typename OUTLETTYPE,
            typename MOMENTTYPE,
            typename TURBULENCETYPE,
            typename ERRORTYPE,
            typename NONDYM,
            typename MOMENTTIMEAVGTYPE>

class SolverTurbulentLES {
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;

    using VectorType = LBMTraits::VectorType;
    using VectorTypeProbeCuda = LBMTraits::VectorTypeProbeCuda;
    using VectorTypeProbeHost = LBMTraits::VectorTypeProbeHost;
    using VectorTypeInt = LBMTraits::VectorType;

    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

public:

    SolverTurbulentLES() = delete;

    SolverTurbulentLES(LBMConstantsPointer &Constants_,
           LBMDataPointer &Data_) {

        Constants = Constants_;
        Data = Data_;

    }

    void initializeSimulation( bool verbose) {

        Constants->Nvel = MODELTYPE::numberOfDiscreteVelocities;

        if (verbose) {
            std::cout << "Number of discrete vel = "<< Constants->Nvel<<".\n";
        }

        initializeSimulationData( verbose );

        setArraySizes();

        NONDYM::nonDimensionalizeInlet(Data, Constants);

        if (verbose) {
            std::cout << "Inlet non-dimensionalized.\n";
        }

        NONDYM::nonDimensionalizeOutlet(Data, Constants);

        if (verbose) {
            std::cout << "Outlet non-dimensionalized.\n";
        }


        INITIALIZATIONTYPE::initialization(Data, Constants);

        if (verbose) {
            std::cout << "Density, Velocity, Distribution functions initialized.\n";
        }

        TURBULENCETYPE::omega(Data, Constants);

        if (verbose) {
            std::cout << "Omega functions initialized.\n";
        }



        std::cout << "Initialization of simulation done." << std::endl;
    }



    void probeWrite(){

        VectorTypeProbeCuda HolderCUDA(4);
        VectorTypeProbeHost HolderHOST(4);

        auto u_view = Data->u.getView();
        auto rho_view = Data->rho.getView();

        int x = Constants -> ProbeLocationLat.x();
        int y = Constants -> ProbeLocationLat.y();
        int z = Constants -> ProbeLocationLat.z();

        auto HoCuView = HolderCUDA.getView();
        auto init = [ = ] __cuda_callable__( int i ) mutable
        {
            HoCuView[0] = u_view(x,y,z).x();
            HoCuView[1] = u_view(x,y,z).y();
            HoCuView[2] = u_view(x,y,z).z();
            HoCuView[3] = rho_view(x,y,z);
        };
        parallelFor< TNL::Devices::Cuda >( 0, 1, init );

        HolderHOST = HolderCUDA;

        auto HoHoView = HolderHOST.getConstView();

        outputerVTK::probeWriteOut(Constants, HoHoView[0], HoHoView[1], HoHoView[2], HoHoView[3]);
    }

    void runSimulation() {

        timer_loop.start();

        int k = 0;
        while(k<Constants -> iterations)
        {

            timer_momentsUpdate.start();
            TURBULENCETYPE::omega(Data, Constants);
            timer_momentsUpdate.stop();

            timer_dumping.start();
            OUTLETTYPE::outletOmega(Data, Constants);
            timer_dumping.stop();


            timer_collision.start();
                COLLISIONTYPE::collision(Data, Constants);
            timer_collision.stop();

            timer_streaming.start();
                STREAMINGTYPE::streaming(Data, Constants);
            timer_streaming.stop();

            timer_bounceback.start();
                BOUNCEBACKWALLTYPE::bounceBackWall(Data, Constants);
                INLETTYPE::inlet(Data, Constants);
                OUTLETTYPE::outlet(Data, Constants);
                SYMMETRYTYPE::symmetry(Data, Constants);
                PERIODICTYPE::periodic(Data, Constants);
            timer_bounceback.stop();

            timer_momentsUpdate.start();
                MOMENTTYPE::momentUpdate(Data, Constants);
            timer_momentsUpdate.stop();

            if(k > (Constants->iterations - Constants->iterationsMomentAvg))
            {
                timer_timeAvg.start();
                    MOMENTTIMEAVGTYPE::momentAdd(Data, Constants);
                timer_timeAvg.stop();
            }

            if(k%Constants->probe_every_it==0 && Constants->probe_every_it > 0 && k>=(Constants->iterations - Constants->probe_iterations))
            {
                timer_output.start();
                probeWrite();
                timer_output.stop();
            }

            if(k%Constants -> err_every_it==0 && k!=0)
            {
                timer_err.start();
                    ERRORTYPE::errorEvaluation(Data, Constants);
                    printf("\n err=%e k=%d\n",Constants->err, k);
                    if (std::isnan(Constants->err))
                    {
                        std::cout << "\n Error is NaN, breaking out.\n";
                        break;
                    }
                timer_err.stop();
            }

            if(k%Constants -> plot_every_it==0 && k )//> (Constants->iterations - Constants->iterationsMomentAvg))
            {

                timer_output.start();
                outputerVTK::variablesVTK(Data, Constants, k/Constants -> plot_every_it, 1);
                timer_output.stop();
            }

            k++;

        }

        //Time averaging
        if(Constants->timeAveraged == false)
        {
            timer_timeAvg.start();
            MOMENTTIMEAVGTYPE::momentAvg(Data, Constants);
            timer_timeAvg.stop();

            timer_output.start();
                outputerVTK::variablesTimeAvgVTK(Data, Constants, 1);
            timer_output.stop();
        }

        timer_loop.stop();
    };

    void convertToLattice(bool verbose) {

        NONDYM::nonDimensionalize(Data, Constants);

    }

    void initializeSimulationData( bool verbose ){



        if( Constants -> plot_every_it  == -1 && Constants -> plot_every != -1) {
            Constants->plot_every_it = std::ceil(Constants->plot_every / Constants->Ct);
        }
        else if( Constants -> plot_every_it  == -1 && Constants -> plot_every == -1) {
            std::cout << "\n !!! CANT GIVE BOTH err_every AND plot_every_it !!!\n";
            assert(false);
        }

        if( Constants -> err_every_it  == -1 && Constants -> err_every != -1) {
            Constants->err_every_it = std::ceil(Constants->err_every / Constants->Ct);
        }
        else if( Constants -> err_every_it  == -1 && Constants -> err_every == -1) {
            std::cout << "\n !!! CANT GIVE BOTH err_every AND err_every_it !!!\n";
            assert(false);
        }
        Constants->iterations = std::ceil(Constants->time / Constants->Ct);

        if (verbose) {
            std::cout << "\nPlotting every " << Constants->plot_every_it << " iterations.\n";
            std::cout << "\nCalculation will run for " << Constants->iterations << " iterations.\n";
        }

        if(Constants->probe_every_it != -1) {
            RealType xProbe = Constants ->ProbeLocation.x();
            RealType yProbe = Constants ->ProbeLocation.y();
            RealType zProbe = Constants ->ProbeLocation.z();

            int xLat = std::round(xProbe/(Constants -> BBmaxx - Constants -> BBminx + 2*(Constants->additional_factor/Constants->resolution_factor))*Constants -> dimX_int);
            int yLat = std::round(yProbe/(Constants -> BBmaxy - Constants -> BBminy+ 2*(Constants->additional_factor/Constants->resolution_factor))*Constants -> dimY_int);
            int zLat = std::round(zProbe/(Constants -> BBmaxz - Constants -> BBminz+ 2*(Constants->additional_factor/Constants->resolution_factor))*Constants -> dimZ_int);


            if(verbose){
                std::cout << "Probe location in lattice " << xLat << ", " << yLat << ", " << zLat << std::endl;
            }

            assert(xLat >= Constants -> dimX_int && "Probe must be in the domain!");
            assert(yLat >= Constants -> dimY_int && "Probe must be in the domain!");
            assert(zLat >= Constants -> dimZ_int && "Probe must be in the domain!");

            VectorTypeInt Probe(xLat, yLat, zLat);
            Constants -> ProbeLocationLat = Probe;

        }
    }

    void setArraySizes(){

        Data->rho.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);
        Data->p.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);
        Data->u.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);

        Data->rhoTimeAvg.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);
        Data->uTimeAvg.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);

        Data->omega.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);

        Data->rho_out.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);
        Data->p_out.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);
        Data->u_out.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);

        Data->u_error.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);

        Data->df.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int, Constants->Nvel);
        Data->df_post.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int, Constants->Nvel);

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
    Timer timer_dumping;
    Timer timer_timeAvg;

};

#endif //SOLVER_H
