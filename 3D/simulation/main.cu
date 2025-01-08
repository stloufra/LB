
#include <iostream>
#include <TNL/Timer.h>
#include <TNL/Logger.h>

#include "src/traits/LBMTraits.h"

#include "src/geometry/geometryMesherBoundary.h"
#include "src/solvers/SolverTurbulentLES.h"

#include "src/postprocesors/outputerVTK.h"
#include "src/postprocesors/outputerMesh.h"
#include "src/postprocesors/outputerStats.h"

#include "src/solvers/models/D3Q27/D3Q27.h"
#include "src/solvers/models/D3Q19/D3Q19.h"
#include "src/solvers/models/D3Q15/D3Q15.h"


using namespace TNL;
using DeviceType = typename LBMTraits::DeviceType;
using VectorType = typename LBMTraits::VectorType;
using RealType = typename LBMTraits::RealType;
using DeviceTypeHost = typename LBMTraits::DeviceTypeHost;
using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

int main() {

    //------------------------INITIALIZATION--------------------------//

    bool runSim = true;

    // initialize data carrier objects
    LBMConstantsPointer Constants;
    LBMDataPointer Data;

    // model types selection
    using Model = D3Q27;

    using Initialisation        = InitializationEquilibriumConstVector<Model>;
    using Collision             = CollisionCumD3Q27Turbulent2015<Model>;
    using Streaming             = StreamingABpush<Model>;
    using BounceBackWall        = BounceBackWallHalf<Model>;
    using Symmetry              = NoSymmetry<Model>;
    using Periodic              = Periodic<Model>;

    using Inlet                 = InletVelocityEquilibrium<Model>;
    using Outlet                = OutletDensityInterpolatedOmegaD3Q27<Model>;
    using Moments               = MomentDensityVelocityN27<Model>;  // SAME AS MODEL NUMBER
    using Error                 = ErrorQuadratic<Model>;
    using Turbulence            = OmegaLES<Model>;
    using NonDim                = NonDimensiolnaliseFactorsVelocity<Model>;
    using TimeAvg               = MomentTimeAvg<Model>;


    //initialize methodical classes
    geometryMesherBoundary Mesher(Constants,
                                  Data);

    SolverTurbulentLES< Model,
            Initialisation,
            Collision,
            Streaming,
            BounceBackWall,
            Symmetry,
            Periodic,
            Inlet,
            Outlet,
            Moments,
            Turbulence,
            Error,
            NonDim,
            TimeAvg> Solver( Constants,
                            Data);



    //------------------------DATA IN--------------------------//


    //set geometry objects -1 streaming from it | no-bounce back - INLET
    //set geometry objects -2 streaming from and into it | bounce back - OUTLET && SYMMETRY
    //set geometry objects -3 periodic

    //resolution 3
    geometryObjectCuboid cuboidInlet({-0.0994f, 0.15f, -0.01f},
                                      {-0.0994f, -0.01f, 0.11f},
                                      {-0.11, 0.15f, -0.01f},-1);


    geometryObjectCuboid cuboidOutlet({0.4595f, 0.15f, -0.01f},
                                      {0.4595f, -0.01f, 0.11f},
                                      {0.46f, 0.15f, -0.01f},-2);

    geometryObjectCuboid cuboidPeriodic1({-0.0994f, 0.0005f, -0.01f},
                                      {-0.0994f, -0.001f, 0.11f},
                                      {0.4595f, 0.0005f, -0.01f},-3);

    geometryObjectCuboid cuboidPeriodic2({-0.0994f, 0.0075f, -0.01f},
                                      {-0.0994f, 0.0075f, 0.11f},
                                      {0.4595f, 0.008f, -0.01f},-3);


    VectorType NormalInlet(-1.f, 0.f, 0.f);

    VectorType velocityInletUniform(72.f, 0.f, 0.f);

    VectorType NormalOutlet(1.f, 0.f, 0.f);

    VectorType NormalPeriodic1(0.f, -1.f, 0.f);

    VectorType NormalPeriodic2(0.f, 1.f, 0.f);

    RealType meanVelocityInlet = 72.f;

    //set simulation initialization
    VectorType Init(60.f, 0.f, 0.f); //change to 1 in z
    Constants->VelocityInit = Init;

    //set physical data
    Constants->conversion_factor_fyz = 1.0f;            //convert to m
    Constants->rho_fyz = 1.293f;                        //[kg/m3]
    Constants->ny_fyz = 2e-5f;                      	//[m2/s]
    Constants->u_guess_fyz = 5.f*meanVelocityInlet;     //[m/s]
    Constants->U_inf = meanVelocityInlet;
    Constants -> omegaDumpingLow = 0.f;
    Constants -> omegaDumpingHigh = 0.004f;		//dumping tau outlet data

    //set lattice data
    Constants->U_lb = 0.09f;                  			// max 0.1

    // set simulation parameters

    Constants->time = 0.1f;                      		//[s]
    Constants->plot_every = 0.1f;               		//[s]
    Constants->err_every = 0.0001f;              		//[s]
    Constants->MomentAvgStart = 0.08f;      		        //[s]
    Constants->MomentAvg_every = 0.02f;      			//[s]

    // set sampling parameters
    VectorType Probe(0.32f, 0.004f, 0.05f);
    Constants->probe_every_it = 1;
    Constants->probe_iterations = 1e4;
    Constants->ProbeLocation = Probe;

    //----------------------LOADING MESH------------------------------//

    outputerMesh::MeshMatrixIn(Data, Constants, "BackwardStepTurbulent", 1);

    //----------------------MESHING GEOMETRY--------------------------//


    Mesher.meshingBoundaryWall(0);

    Mesher.meshingBoundaryConditionPeriodic(cuboidPeriodic1,NormalPeriodic1, 16, 1);
    Mesher.meshingBoundaryConditionPeriodic(cuboidPeriodic2,NormalPeriodic2, 1, 1);
    Mesher.meshingBoundaryConditionInletUniform(cuboidInlet, NormalInlet, velocityInletUniform,0);
    Mesher.meshingBoundaryConditionOutlet(cuboidOutlet, NormalOutlet, Constants->rho_fyz,1);

    Mesher.arrayTransfer(1);



    //----------------------MESHING OUTPUT--------------------------//


    outputerVTK::MeshVTK(Data, Constants, "meshIN");



    //----------------------SOLVING PROBLEM------------------------//


    Solver.convertToLattice(1);
    Solver.initializeSimulation(1);

    if(runSim) {


        Solver.runSimulation();
    }

    //----------------------TIMERS OUTPUT--------------------------//

	outputerStats::Statistics(Solver, Constants, 0, 1);

    return 0;
}
