#pragma hd_warning_disable
#pragma nv_exec_check_disable


#include <iostream>
#include <TNL/Timer.h>
#include <TNL/Logger.h>

#include "src/geometry/geometryMesherBoundary.h"
#include "src/geometry/geometryObjectCuboid.h"
#include "src/solvers/Solver.h"
#include "src/traits/LBMTraits.h"
#include "src/postprocesors/outputerVTK.h"
#include "src/postprocesors/outputerMesh.h"

#include "src/solvers/models/D3Q27/D3Q27.h"
#include "src/solvers/models/D3Q19/D3Q19.h"
#include "src/solvers/models/D3Q15/D3Q15.h"


using namespace TNL;
using DeviceType = typename LBMTraits::DeviceType;
using VectorType = typename LBMTraits::VectorType;
using DeviceTypeHost = typename LBMTraits::DeviceTypeHost;
using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

int main() {


    //------------------------INITIALIZATION--------------------------//


    // initialize data carrier objects
    LBMConstantsPointer Constants;
    LBMDataPointer Data;

    // model types selection
    using Model = D3Q15;

    using Initialisation = InitializationEquilibriumVariables<Model>;
    using Collision = CollisionSRT<Model>;
    using Streaming = StreamingAB<Model>;
    using BounceBackWall = BounceBackWallHalf<Model>;
    using Inlet = InletVelocity<Model>;
    using Outlet = OutletDensityEquilibrium<Model>;
    using Moments  = MomentDensityVelocityN15<Model>;  // SAME AS MODEL NUMBER
    using Error = ErrorQuadratic<Model>;
    using NonDim = NonDimensiolnaliseFactorsVelocity<Model>;


    //initialize timers
    Timer timerMeshingBoundary;
    Logger logger(50, std::cout);

    //initialize methodical classes
    geometryMesherBoundary Mesher(Constants,
                                  Data);

    Solver< Model,
            Initialisation,
            Collision,
            Streaming,
            BounceBackWall,
            Inlet,
            Outlet,
            Moments,
            Error,
            NonDim> Solver( Constants,
                            Data);


    //------------------------DATA IN--------------------------//

    //set simulation initialization
    VectorType Init(1.f, 0.f, 0.f);
    Constants->VelocityInit = Init;
    Constants->InitFileName = "variablesLattice199-backup-2-2-factor";


    //set meshing data
    Constants->resolution_factor = 3.f;                              // needs to be 1 or greater integer
    Constants->additional_factor = 2.f;                              // at least 1 for additional wall around
    Constants->point_outside = {0.f, 0.f, 20.f};
    Constants->file_name = "Dummy.off";

    //set geometry objects

    //resolution 3
    geometryObjectCuboid cuboidInlet1({15.f, 160.f, -80.f},
                                      {-15.f, 160.f, -80.f},
                                      {15.f, 120.f, -79.5f},
                                      3);

    geometryObjectCuboid cuboidInlet2({15.f, 200.f, 15.f},
                                      {-15.f, 200.f, 15.f},
                                      {15.f, 199.5f, -5.f},
                                      4);

    geometryObjectCuboid cuboidOutlet({15.f, 0.f, 15.f},
                                      {-15.f, 0.f, 15.f},
                                      {15.f, 0.3f, -15.f},
                                      5);


    VectorType VelocityInlet1(0.f, 0.f, 0.1f);
    VectorType VelocityInlet2(0.f, -0.2f, 0.f);
    VectorType NormalInlet1(0.f, 0.f, -1.f);
    VectorType NormalInlet2(0.f, 1.f, 0.f);
    VectorType NormalOutlet(0.f, -1.f, 0.f);


    //set physical data
    Constants->rho_fyz = 1000.f;                      //[kg/m3]
    Constants->ny_fyz = 10e-5f;                       //[m2/s]
    Constants->u_guess_fyz = 0.5f;                   //[m/s] //TODO should be automatically calculated
    Constants->Fx_fyz = 10.f;                         //[kg/m3/s2]  <- force density
    Constants->Fy_fyz = 0.0f;                         //[kg/m3/s2]  <- force density
    Constants->Fz_fyz = 0.0f;                         //[kg/m3/s2]  <- force density
    Constants->conversion_factor_fyz = 1.0f / 1000.f;    // convert to m

    //set lattice data

    Constants->U_lb = 0.09f;                  // max 0.1 (Book suggests max 0.2)

    // set simulation parameters

    Constants->time = 2.f;               //[s]
    Constants->plot_every = 0.01f;         //[s]
    Constants->err_every = 0.002f;         //[s]

    //----------------------LOADING MESH------------------------------//

    outputerMesh::MeshMatrixIn(Data, Constants, "mesh", 1);

    //----------------------MESHING GEOMETRY--------------------------//

    timerMeshingBoundary.start();
    Mesher.meshingBoundaryWall(0);
    Mesher.meshingBoundaryConditionInlet(cuboidInlet1, NormalInlet1, VelocityInlet1, 1);
    Mesher.meshingBoundaryConditionInlet(cuboidInlet2, NormalInlet2, VelocityInlet2, 1);
    Mesher.meshingBoundaryConditionOutlet(cuboidOutlet, NormalOutlet, 1000.f,
                                          1); //if density - 1 then density is from noditself
    Mesher.compileBoundaryArrayInlets(1);
    Mesher.compileBoundaryArrayOutlets(1);
    Mesher.arrayTransfer(1);
    timerMeshingBoundary.stop();


    //----------------------MESHING OUTPUT--------------------------//


    outputerVTK::MeshVTK(Data, Constants, "meshIN");

    //----------------------SOLVING PROBLEM------------------------//

    Solver.convertToLattice(1);
    Solver.initializeSimulation(1);
    Solver.runSimulation();

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


    return 0;
}


