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
#include "src/solvers/models/D3Q27.h"

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

    //TODO: using Initialization = initializationEqulibrium;
    using ModelData = D3Q27Data;
    using Model = D3Q27<ModelData>; //TODO <Initialization >;
    //initialize timers
    Timer timerMeshingBoundary;
    Logger logger(50, std::cout);

    //initialize methodical classes


    geometryMesherBoundary Mesher(Constants,
                                  Data);

    Solver<Model, ModelData> Solver(Constants,
                  Data);
    //------------------------DATA IN--------------------------//

    //set meshing data
    Constants->resolution_factor = 3;                              // needs to be 1 or greater integer
    Constants->additional_factor = 2;                              // at least 1 for additional wall around
    Constants->point_outside = {0, 0, 20};
    Constants->file_name = "Dummy.off";

    //set geometry objects

    //resolution 3
    geometryObjectCuboid cuboidInlet1({15, 160, -80},
                                      {-15, 160, -80},
                                      {15, 120, -79.5},
                                      3);

    geometryObjectCuboid cuboidInlet2({15, 200, 15},
                                      {-15, 200, 15},
                                      {15, 199.5, -5},
                                      4);

    geometryObjectCuboid cuboidOutlet({15, 0, 15},
                                      {-15, 0, 15},
                                      {15, 0.3, -15},
                                      5);



    VectorType VelocityInlet1( 0, 0, 0.1);
    VectorType VelocityInlet2(0, -0.2, 0);
    VectorType NormalInlet1(0, 0, -1);
    VectorType NormalInlet2(0, 1, 0);
    VectorType NormalOutlet(0, -1, 0);

    VectorType VelocityInit( 0, 0, 0);

    //set physical data
    Constants->rho_fyz = 1000.f;                      //[kg/m3]
    Constants->ny_fyz = 10e-5f;                       //[m2/s]
    Constants->u_guess_fyz = 0.5f;                   //[m/s]
    Constants->Fx_fyz = 10.f;                         //[kg/m3/s2]  <- force density
    Constants->Fy_fyz = 0.0f;                         //[kg/m3/s2]  <- force density
    Constants->Fz_fyz = 0.0f;                         //[kg/m3/s2]  <- force density
    Constants->conversion_factor_fyz = 1.0 / 1000;    // convert to m

    //set lattice data

    Constants->U_lb = 0.09;                  // max 0.1 (Book suggests max 0.2)

    // set simulation parameters

    Constants -> time =0.01f;               //[s]
    Constants -> plot_every=0.001f;         //[s]

    //----------------------LOADING MESH------------------------------//

    outputerMesh::MeshMatrixIn(Data, Constants, "mesh", 1);

    //----------------------MESHING GEOMETRY--------------------------//

    timerMeshingBoundary.start();
        Mesher.meshingBoundaryWall(0);
        Mesher.meshingBoundaryConditionInlet(cuboidInlet1, NormalInlet1, VelocityInlet1, 1);
        Mesher.meshingBoundaryConditionInlet(cuboidInlet2, NormalInlet2, VelocityInlet2, 1);
        Mesher.meshingBoundaryConditionOutlet(cuboidOutlet, NormalOutlet, 0.788, 1); //if density - 1 then density is from noditself
        Mesher.compileBoundaryArrayInlets(1);
        Mesher.compileBoundaryArrayOutlets(1);
        Mesher.arrayTransfer(1);
    timerMeshingBoundary.stop();


    //----------------------MESHING OUTPUT--------------------------//


        outputerVTK::MeshVTK(Data, Constants, "meshIN");


    //----------------------SOLVING PROBLEM------------------------//

        Solver.convertToLattice(1);
        Solver.initializeSimulation(VelocityInit, 1);
        //outputerVTK::variablesLatticeVTK(Data, Constants, -1, 1);
        Solver.runSimulation();
        outputerVTK::distributionFunctionVTK(Data, Constants, 3, 1);

    //----------------------TIMERS OUTPUT--------------------------//


    logger.writeHeader("Timing of sections 1) Whole loop");
    //logger.writeSystemInformation(true);
    Solver.timer_loop.writeLog( logger, 0 );
    logger.writeSeparator();
    logger.writeHeader("Collision");
    Solver.timer_collision.writeLog( logger, 0 );
    logger.writeSeparator();
    logger.writeHeader("Streaming");
    Solver.timer_streaming.writeLog( logger, 0 );
    logger.writeSeparator();
    logger.writeHeader("Bounce back");
    Solver.timer_bounceback.writeLog( logger, 0 );
    logger.writeSeparator();
    logger.writeHeader("Moments Update");
    Solver.timer_momentsUpdate.writeLog( logger, 0 );
    logger.writeSeparator();
    logger.writeHeader("Error Calculation");
    Solver.timer_err.writeLog( logger, 0 );
    logger.writeSeparator();
    logger.writeHeader("Writting Output");
    Solver.timer_output.writeLog( logger, 0 );


    return 0;
}


