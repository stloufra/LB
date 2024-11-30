#ifndef GEOMETRYMESHEROFF_H
#define GEOMETRYMESHEROFF_H

#include <string>
#include <omp.h>
#include <vector>


#include "geometryHandlerOFF.h"
#include "geometryObjectCuboid.h"
#include "../traits/LBMTraits.h"
#include "../data/LBMData.h"
//#include "../data/LBMConstants.h"

#include <TNL/Algorithms/parallelFor.h>

#pragma once

using namespace TNL;
using namespace TNL::Algorithms;

class geometryMesherOFF {
public:

    using RealType = typename LBMTraits::RealType;
    using DeviceType = typename LBMTraits::DeviceType;
    using DeviceTypeHost = typename LBMTraits::DeviceTypeHost;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;
    using ArrayTypeMeshHost = typename LBMTraits::ArrayTypeMeshHost;
    using ArrayTypeBoundaryInletHost = typename LBMTraits::ArrayTypeBoundaryInletHost;
    using ArrayTypeBoundaryOutletHost = typename LBMTraits::ArrayTypeBoundaryOutletHost;
    using geometryHandlerOFFPointer = TNL::Pointers::SharedPointer<geometryHandlerOFF, DeviceTypeHost>;


    geometryMesherOFF() = delete;

    geometryMesherOFF(LBMConstantsPointer Constants_,
                      LBMDataPointer Data_) {

        Constants = Constants_;
        Data = Data_;
    }

    void initialization(bool verbose) {

        Constants->dimX_int = std::ceil(
                (Constants->BBmaxx - Constants->BBminx) * Constants->resolution_factor +
                2 * Constants->additional_factor); // *-|-*--*--*... inside same nm. of points and spaces
        Constants->dimY_int = std::ceil(
                (Constants->BBmaxy - Constants->BBminy) * Constants->resolution_factor +
                2 * Constants->additional_factor);
        Constants->dimZ_int = std::ceil(
                (Constants->BBmaxz - Constants->BBminz) * Constants->resolution_factor +
                2 * Constants->additional_factor);

        Data->meshFluidHost.setSizes(Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int);

        auto mesh_view = Data->meshFluidHost.getView();

        auto init_mesh = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i  ) mutable
        {
            mesh_view(i.x(), i.y(), i.z()) = 0;
        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceTypeHost>(begin, end, init_mesh);

        if (verbose) {
            std::cout << "Resolution factor: " << Constants->resolution_factor << std::endl;
            std::cout << "Additional factor: " << Constants->additional_factor << std::endl;
            std::cout << "Dimension X: " << Constants->dimX_int << std::endl;
            std::cout << "Dimension Y: " << Constants->dimY_int << std::endl;
            std::cout << "Dimension Z: " << Constants->dimZ_int << std::endl;
        }
    }

    void meshingPolyhedron(geometryHandlerOFFPointer &geometry_handler, d3 pnt_outside, bool verbose=1) {

        RealType lookx, looky, lookz;

        auto mesh_view = Data->meshFluidHost.getView();

        // res.f. = 10 | BminX-BmaxX=10 | add.f = 10
        // -0.05    0.05  0.15  9.95   10.05
        //   * - | - * - - * ... * - | - *

        const int totalIterations = Constants->dimX_int * Constants->dimY_int * Constants->dimZ_int;
        int progressCounter = 0;

        if (verbose)
        {
            std::cout << "\nMeshing Polyhedron:" <<std::endl;
        }

        for (int i = 0; i < Constants->dimX_int; i++) {
            for (int j = 0; j < Constants->dimY_int; j++) {
                for (int k = 0; k < Constants->dimZ_int; k++) {
                    lookx = (static_cast<double>(i)+0.5 - Constants->additional_factor) / Constants->resolution_factor + Constants->BBminx; // 0.5 for halfway between wall and first node
                    looky = (static_cast<double>(j)+0.5 - Constants->additional_factor) / Constants->resolution_factor + Constants->BBminy;
                    lookz = (static_cast<double>(k)+0.5 - Constants->additional_factor) / Constants->resolution_factor + Constants->BBminz;


                    d3 pnt_ask = {lookx, looky, lookz};

                    if (geometry_handler->polyhedronInside(pnt_ask, pnt_outside)) {
                        mesh_view(i, j, k) = 1;
                    }

                   progressCounter++;
                }
            }
            if (verbose)
            {
                if (progressCounter % (totalIterations / 100) <= (totalIterations / 1000)|| progressCounter == totalIterations)
                {
                    printProgressBar(progressCounter, totalIterations);
                }
            }
        }

        if (verbose)
        {
            std::cout <<"\n" << std::endl;
        }

        addOutsideLayer();


        if (verbose) {
            printBoundaryPositions();
            printCuboidsForBC();

            std::cout << "Number of lattice points is "
                      << Constants->dimX_int * Constants->dimY_int * Constants->dimZ_int << "\n";
        }
    }

    void printProgressBar(int current, int total, int barWidth = 50){
        float progress = static_cast<float>(current) / total;
        int position = static_cast<int>(barWidth * progress);

        std::cout << "\r|";
        for (int i = 0; i < barWidth; ++i) {
            if (i < position)
                std::cout << "O";
            else
                std::cout << "-";
        }
        std::cout << "| " << static_cast<int>(progress * 100) << "%";
        std::cout.flush();
    }

    void printBoundaryPositions() {
        auto calculatePosition = [&](int index, double minBound, double additionalFactor, double resolutionFactor) {
            return (static_cast<double>(index) + 0.5 - additionalFactor) / resolutionFactor + minBound;
        };

        std::cout << "Boundary positions:\n";

        // X-direction
        std::cout << "X-direction:\n";
        for (int i : {0, 1, 2, Constants->dimX_int - 3, Constants->dimX_int  - 2,Constants->dimX_int - 1}) {
            double x = calculatePosition(i, Constants->BBminx, Constants->additional_factor, Constants->resolution_factor);
            std::cout << "  Index " << i << " -> Position: " << x << "\n";
        }

        // Y-direction
        std::cout << "Y-direction:\n";
        for (int j : {0, 1, 2, Constants->dimY_int  - 3, Constants->dimY_int - 2, Constants->dimY_int - 1}) {
            double y = calculatePosition(j, Constants->BBminy, Constants->additional_factor, Constants->resolution_factor);
            std::cout << "  Index " << j << " -> Position: " << y << "\n";
        }

        // Z-direction
        std::cout << "Z-direction:\n";
        for (int k : {0, 1, 2, Constants->dimZ_int - 3, Constants->dimZ_int - 2, Constants->dimZ_int - 1}) {
            double z = calculatePosition(k, Constants->BBminz, Constants->additional_factor, Constants->resolution_factor);
            std::cout << "  Index " << k << " -> Position: " << z << "\n";
        }
    }

    void addOutsideLayer(){

        auto mesh_view = Data->meshFluidHost.getView();


        if(Constants -> additional_factor >0)
        {

            // Set boundary points to 0
            for (int i = 0; i < Constants->dimX_int; i++) {
                for (int j = 0; j < Constants->dimY_int; j++) {
                    // Set boundaries along Z axis
                    mesh_view(i, j, 0) = 0;
                    mesh_view(i, j, Constants->dimZ_int - 1) = 0;
                }
            }

            for (int i = 0; i < Constants->dimX_int; i++) {
                for (int k = 0; k < Constants->dimZ_int; k++) {
                    // Set boundaries along Y axis
                    mesh_view(i, 0, k) = 0;
                    mesh_view(i, Constants->dimY_int - 1, k) = 0;
                }
            }

            for (int j = 0; j < Constants->dimY_int; j++) {
                for (int k = 0; k < Constants->dimZ_int; k++) {
                    // Set boundaries along X axis
                    mesh_view(0, j, k) = 0;
                    mesh_view(Constants->dimX_int - 1, j, k) = 0;
                }
            }
        }
    }


    void printCuboidsForBC() {

    auto calculatePosition = [&](int index, double minBound, double additionalFactor, double resolutionFactor) {
        return (static_cast<double>(index) + 0.5 - additionalFactor) / resolutionFactor + minBound;
    };

    auto printCuboid = [&](const std::string& tag, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max) {
        std::cout << std::fixed << std::setprecision(6)
                  << "geometryObjectCuboid " << tag << "({"
                  << x_min << "f, " << y_min << "f, " << z_max << "f},\n"
                  << "                                      {"
                  << x_min << "f, " << y_max << "f, " << z_max << "f},\n"
                  << "                                      {"
                  << x_max << "f, " << y_min << "f, " << z_min << "f}, \"" << tag << "\");\n";
    };

    // Max values
    double x_min = calculatePosition(0, Constants->BBminx, Constants->additional_factor, Constants->resolution_factor);
    double x_max = calculatePosition(Constants-> dimX_int- 1, Constants->BBminx, Constants->additional_factor, Constants->resolution_factor);
    double y_min = calculatePosition(0, Constants->BBminy, Constants->additional_factor, Constants->resolution_factor);
    double y_max = calculatePosition(Constants-> dimY_int- 1, Constants->BBminy, Constants->additional_factor, Constants->resolution_factor);
    double z_min = calculatePosition(0, Constants->BBminz, Constants->additional_factor, Constants->resolution_factor);
    double z_max = calculatePosition(Constants-> dimZ_int - 1, Constants->BBminz, Constants->additional_factor, Constants->resolution_factor);


    // X-Direction Inlet and Outlet
    double x0 = calculatePosition(0, Constants->BBminx, Constants->additional_factor, Constants->resolution_factor);
    double x1 = calculatePosition(1, Constants->BBminx, Constants->additional_factor, Constants->resolution_factor);
    double x2 = calculatePosition(2, Constants->BBminx, Constants->additional_factor, Constants->resolution_factor);
    double x_in1 = (x0 + x1) / 2.0;
    double x_in2 = (x1 + x2) / 2.0;

    double x3 = calculatePosition(Constants-> dimX_int- 3, Constants->BBminx, Constants->additional_factor, Constants->resolution_factor);
    double x4 = calculatePosition(Constants-> dimX_int- 2, Constants->BBminx, Constants->additional_factor, Constants->resolution_factor);
    double x5 = calculatePosition(Constants-> dimX_int- 1, Constants->BBminx, Constants->additional_factor, Constants->resolution_factor);
    double x_out1 = (x3 + x4) / 2.0;
    double x_out2 = (x4 + x5) / 2.0;

    printCuboid("cuboidX_in", x_in1, x_in2, y_min, y_max, z_min, z_max);
    printCuboid("cuboidX_out", x_out1, x_out2, y_min, y_max, z_min, z_max);

    // Y-Direction Inlet and Outlet
    double y0 = calculatePosition(0, Constants->BBminy, Constants->additional_factor, Constants->resolution_factor);
    double y1 = calculatePosition(1, Constants->BBminy, Constants->additional_factor, Constants->resolution_factor);
    double y2 = calculatePosition(2, Constants->BBminy, Constants->additional_factor, Constants->resolution_factor);
    double y_in1 = (y0 + y1) / 2.0;
    double y_in2 = (y1 + y2) / 2.0;

    double y3 = calculatePosition(Constants-> dimY_int- 3, Constants->BBminy, Constants->additional_factor, Constants->resolution_factor);
    double y4 = calculatePosition(Constants-> dimY_int- 2, Constants->BBminy, Constants->additional_factor, Constants->resolution_factor);
    double y5 = calculatePosition(Constants-> dimY_int- 1, Constants->BBminy, Constants->additional_factor, Constants->resolution_factor);
    double y_out1 = (y3 + y4) / 2.0;
    double y_out2 = (y4 + y5) / 2.0;

    printCuboid("cuboidY_in", x0, x2, y_in1, y_in2, z_min, z_max);
    printCuboid("cuboidY_out", x0, x2, y_out1, y_out2, z_min, z_max);

    // Z-Direction Inlet and Outlet
    double z0 = calculatePosition(0, Constants->BBminz, Constants->additional_factor, Constants->resolution_factor);
    double z1 = calculatePosition(1, Constants->BBminz, Constants->additional_factor, Constants->resolution_factor);
    double z2 = calculatePosition(2, Constants->BBminz, Constants->additional_factor, Constants->resolution_factor);
    double z_in1 = (z0 + z1) / 2.0;
    double z_in2 = (z1 + z2) / 2.0;

    double z3 = calculatePosition(Constants-> dimZ_int- 3, Constants->BBminz, Constants->additional_factor, Constants->resolution_factor);
    double z4 = calculatePosition(Constants-> dimZ_int- 2, Constants->BBminz, Constants->additional_factor, Constants->resolution_factor);
    double z5 = calculatePosition(Constants-> dimZ_int- 1, Constants->BBminz, Constants->additional_factor, Constants->resolution_factor);
    double z_out1 = (z3 + z4) / 2.0;
    double z_out2 = (z4 + z5) / 2.0;

    printCuboid("cuboidZ_in", x0, x2, y_min, y_max, z_in1, z_in2);
    printCuboid("cuboidZ_out", x0, x2, y_min, y_max, z_out1, z_out2);
}



    std::vector <boundaryConditionOutlet> boundary_vector_outlet_individual;
    std::vector <boundaryConditionInlet> boundary_vector_inlet_individual;
    std::vector <boundaryConditionInlet> boundary_vector_inlet;
    std::vector <boundaryConditionOutlet> boundary_vector_outlet;

    boundaryConditionInlet boundary_condition_inlet;

    LBMConstantsPointer Constants;
    LBMDataPointer Data;
};


#endif //GEOMETRYMESHER_H

#pragma clang diagnostic pop