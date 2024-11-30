#include <iostream>
#include "src/geometry/geometryHandlerOFF.h"
#include "src/geometry/geometryMesherOFF.h"
#include "src/geometry/geometryObjectCuboid.h"
#include "src/traits/LBMTraits.h"
#include "src/postprocesors/outputerVTK.h"
#include "src/postprocesors/outputerMesh.h"
#include <TNL/Timer.h>
#include <TNL/Logger.h>


using namespace TNL;
using DeviceType = typename LBMTraits::DeviceType;
using VectorType = typename LBMTraits::VectorType;
using DeviceTypeHost = typename LBMTraits::DeviceTypeHost;
using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;
using geometryHandlerOFFPointer = TNL::Pointers::SharedPointer<geometryHandlerOFF, DeviceTypeHost>;

int main() {


    //------------------------INITIALIZATION--------------------------//


    // initialize data carrier objects
    LBMConstantsPointer Constants;
    LBMDataPointer Data;
    geometryHandlerOFFPointer Handler;

    //initialize timers
    Timer timer_handling;              // timer for handling the polyhedron
    Timer timer_meshing;               // timer for meshing the polyhedron
    Timer timer_VTK;                   // timer for writing the VTK file
    Logger logger(50, std::cout);

    geometryMesherOFF Mesher(Constants,
                          Data);

    //set meshing data
    Constants->resolution_factor = 500;
    Constants->additional_factor = 1;                      // at least 1 for additional wall around
    Constants->point_outside = {-1, -1, -1};
    Constants->file_name = "./CAD/TurbulentProfilePipe.off";
    Constants->mesh_name = "PipeProfileFine";



    //----------------------HANDLING GEOMETRY--------------------------//


    timer_handling.start();
        Handler->polyhedronFromFile(Constants, 0);    //read the file and store the polyhedron in the handler object
        Handler->polyhedronBbox(1);                     //compute the bounding box of the polyhedron
        Handler->writeToConstants(Constants);               //write the BBOX into the Constats
    timer_handling.stop();



    //----------------------MESHING GEOMETRY--------------------------//


    timer_meshing.start();
        Mesher.initialization(1);                                       //initialize    the mesh
        Mesher.meshingPolyhedron(Handler, Constants->point_outside, 1);
    timer_meshing.stop();


    //----------------------MESHING OUTPUT--------------------------//

    timer_VTK.start();
        outputerVTK::MeshVTK(Data, Constants, Constants->mesh_name);
    timer_VTK.stop();

    outputerMesh::MeshMatrixOut(Data, Constants, Constants->mesh_name);

    //----------------------TIMERS OUTPUT--------------------------//
    logger.writeHeader("Handling Geometry");
        timer_handling.writeLog(logger, 0);
    logger.writeSeparator();
    logger.writeHeader("Meshing Polyhedron");
        timer_meshing.writeLog(logger, 0);
    logger.writeSeparator();
    logger.writeHeader("Writting Geometry");
        timer_VTK.writeLog(logger, 0);
    logger.writeSeparator();

    printf("Number of lattice points - %d. \n", Constants->dimX_int*Constants->dimY_int*Constants->dimZ_int);

    return 0;
}


