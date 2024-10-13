#pragma hd_warning_disable
#pragma nv_exec_check_disable


#include <iostream>
#include <TNL/Timer.h>
#include <TNL/Logger.h>

#include "src/geometry/geometryMesherBoundary.h"
#include "src/geometry/geometryObjectCuboid.h"
#include "src/solvers/SolverTurbulentLES.h"
#include "src/traits/LBMTraits.h"
#include "src/postprocesors/outputerVTK.h"
#include "src/postprocesors/outputerMesh.h"

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
    using Streaming             = StreamingAB<Model>;
    using BounceBackWall        = BounceBackWallHalf<Model>;
    using Symmetry              = NoSymmetry<Model>;

    using Inlet                 = InletVelocityEquilibrium<Model>;
    using Outlet                = OutletDensityInterpolated<Model>;
    using Moments               = MomentDensityVelocityN27<Model>;  // SAME AS MODEL NUMBER
    using Error                 = ErrorQuadratic<Model>;
    using Turbulence            = OmegaLES<Model>;
    using NonDim                = NonDimensiolnaliseFactorsVelocity<Model>;
    using TimeAvg               = MomentTimeAvg<Model>;


    //initialize timers
    Timer timerMeshingBoundary;
    Logger logger(50, std::cout);

    //initialize methodical classes
    geometryMesherBoundary Mesher(Constants,
                                  Data);

    SolverTurbulentLES< Model,
            Initialisation,
            Collision,
            Streaming,
            BounceBackWall,
            Symmetry,
            Inlet,
            Outlet,
            Moments,
            Turbulence,
            Error,
            NonDim,
            TimeAvg> Solver( Constants,
                            Data);



    //------------------------DATA IN--------------------------//

    //set simulation initialization
    VectorType Init(0.f, 0.f, 3.f); //change to 1 in z
    Constants->VelocityInit = Init;


   /* //set meshing data
    Constants->resolution_factor = 0.1;
    Constants->additional_factor = 6;                              // at least 1 for additional wall around
    Constants->point_outside = {2, 1000, 0};
    Constants->file_name = "Fany.off";*/


    //set geometry objects -1 streaming from it | no-bounce back - INLET
    //set geometry objects -2 streaming from and into it | bounce back - OUTLET

    //resolution 3
    geometryObjectCuboid cuboidInlet({-0.099f, 0.15f, -0.01f},
                                      {-0.099f, -0.01f, 0.11f},
                                      {-0.11, 0.15f, -0.01f},-1);


    geometryObjectCuboid cuboidOutlet({0.4575f, 0.15f, -0.01f},
                                      {0.4575f, -0.01f, 0.11f},
                                      {0.46f, 0.15f, -0.01f},-2);


    VectorType NormalInlet(-1.f, 0.f, 0.f);

    VectorType velocityInletUniform(46.78f, 0.f, 0.f);

    VectorType NormalOutlet(1.f, 0.f, 0.f);

    //inlet parabolic data

    d3 inletCenter = {-0.1f , 0.075, 0.06};
    RealType inletDimX = 0.f;
    RealType inletDimY = 0.15f;
    RealType inletDimZ = 0.1f;
    RealType meanVelocityInlet = 46.78f;       // 5

    //dumping tau outlet data
    Constants -> omegaDumpingLow = 0.f;
    Constants -> omegaDumpingHigh = 0.05f; // tau(Re=1000) = 0.55 -> 0.2


    //set physical data
    Constants->rho_fyz = 1.293f;                        //[kg/m3]     1000
    Constants->ny_fyz = 2*10e-5f;                       //[m2/s]
    Constants->u_guess_fyz = 4.f*meanVelocityInlet;     //[m/s]
    Constants->Fx_fyz = 0.0f;                           //[kg/m3/s2]  <- force density
    Constants->Fy_fyz = 0.0f;                           //[kg/m3/s2]  <- force density
    Constants->Fz_fyz = 0.0f;                           //[kg/m3/s2]  <- force density
    Constants->conversion_factor_fyz = 1.0f;            // convert to m

    //set lattice data

    Constants->U_lb = 0.09f;                  // max 0.1 (Book suggests max 0.2)

    // set simulation parameters

    Constants->time = 80.0f;                      //[s]
    Constants->plot_every = 0.5f;               //[s]
    Constants->err_every = 0.001f;              //[s]
    Constants->iterationsMomentAvg = 10000;      //[1]

    // set sampling parameters
    Constants->probe_every_it = 100;
    VectorType Probe(0.301f, 0.075f, 0.05f);
    Constants->ProbeLocation = Probe;

    //----------------------LOADING MESH------------------------------//

    outputerMesh::MeshMatrixIn(Data, Constants, "BackwardStepTurbulent", 1);

    //----------------------MESHING GEOMETRY--------------------------//

    timerMeshingBoundary.start();
        Mesher.meshingBoundaryWall(0);
        //Mesher.meshingBoundaryConditionInletParaboloidRectangle( cuboidInlet, inletCenter, inletDimX, inletDimY, inletDimZ, NormalInlet, meanVelocityInlet, 1 );


        Mesher.meshingBoundaryConditionInletUniform(cuboidInlet, NormalInlet, velocityInletUniform,0);
        Mesher.meshingBoundaryConditionOutlet(cuboidOutlet, NormalOutlet, Constants->rho_fyz,
                                                1); //TODO: if density = -1 then density is from nod itself
        Mesher.compileBoundaryArrayWall(1);
        Mesher.compileBoundaryArrayInlets(1);
        Mesher.compileBoundaryArrayOutlets(1);
        Mesher.arrayTransfer(1);
    timerMeshingBoundary.stop();


    //----------------------MESHING OUTPUT--------------------------//


    outputerVTK::MeshVTK(Data, Constants, "meshIN");



    //----------------------SOLVING PROBLEM------------------------//


    Solver.convertToLattice(1);
    Solver.initializeSimulation(1);

    if(runSim) {
        Solver.runSimulation();
    }

    //----------------------TIMERS OUTPUT--------------------------//


    logger.writeHeader("Timing of sections 1) Whole loop");
    //logger.writeSystemInformation(true);
    Solver.timer_loop.writeLog(logger, 0);
    logger.writeSeparator();
    logger.writeHeader("Collision");
    Solver.timer_collision.writeLog(logger, 0);
    logger.writeSeparator();
    logger.writeHeader("Streaming");
    Solver.timer_streaming.writeLog(logger, 0);
    logger.writeSeparator();
    logger.writeHeader("Bounce back");
    Solver.timer_bounceback.writeLog(logger, 0);
    logger.writeSeparator();
    logger.writeHeader("Moments Update");
    Solver.timer_momentsUpdate.writeLog(logger, 0);
    logger.writeSeparator();
    logger.writeHeader("Error Calculation");
    Solver.timer_err.writeLog(logger, 0);
    logger.writeSeparator();
    logger.writeHeader("Writting Output");
    Solver.timer_output.writeLog(logger, 0);
    logger.writeSeparator();
    logger.writeHeader("TimeDumping");
    Solver.timer_dumping.writeLog(logger, 0);
    logger.writeSeparator();
    logger.writeHeader("Time Averaging");
    Solver.timer_timeAvg.writeLog(logger, 0);


    return 0;
}
