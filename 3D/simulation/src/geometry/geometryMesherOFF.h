#ifndef GEOMETRYMESHEROFF_H
#define GEOMETRYMESHEROFF_H

#include <string>
#include <omp.h>
#include <vector>


#include "geometryHandlerOFF.h"
#include "geometryMesherBoundary.h"
#include "geometryObjectCuboid.h"
#include "../traits/LBMTraits.h"
//#include "../data/LBMData.h"
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
                2 * Constants->additional_factor * Constants->resolution_factor);
        Constants->dimY_int = std::ceil(
                (Constants->BBmaxy - Constants->BBminy) * Constants->resolution_factor +
                2 * Constants->additional_factor * Constants->resolution_factor);
        Constants->dimZ_int = std::ceil(
                (Constants->BBmaxz - Constants->BBminz) * Constants->resolution_factor +
                2 * Constants->additional_factor * Constants->resolution_factor);

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

    void meshingPolyhedron(geometryHandlerOFFPointer &geometry_handler, d3 pnt_outside, bool verbose) {

        RealType lookx, looky, lookz;

        auto mesh_view = Data->meshFluidHost.getView();


        for (int i = 0; i < Constants->dimX_int; i++) {
            for (int j = 0; j < Constants->dimY_int; j++) {
                for (int k = 0; k < Constants->dimZ_int; k++) {
                    lookx = static_cast<double>(i) / Constants->resolution_factor + Constants->BBminx -
                            Constants->additional_factor;
                    looky = static_cast<double>(j) / Constants->resolution_factor + Constants->BBminy -
                            Constants->additional_factor;
                    lookz = static_cast<double>(k) / Constants->resolution_factor + Constants->BBminz -
                            Constants->additional_factor;

                    d3 pnt_ask = {lookx, looky, lookz};

                    if (geometry_handler -> polyhedronInside(pnt_ask, pnt_outside)) {
                        mesh_view(i, j, k) = 1;
                    }
                }
            }
        }

        /*auto meshPoly = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i  ) mutable
        {

            lookx = static_cast<double>(i.x()) / Constants->resolution_factor + Constants->BBminx -
                    Constants->additional_factor;
            looky = static_cast<double>(i.y()) / Constants->resolution_factor + Constants->BBminy -
                    Constants->additional_factor;
            lookz = static_cast<double>(i.z()) / Constants->resolution_factor + Constants->BBminz -
                    Constants->additional_factor;

            d3 pnt_ask = {lookx, looky, lookz};

            if (geometry_handler -> polyhedronInside(pnt_ask, pnt_outside)) {
                mesh_view(i.x(), i.y(), i.z()) = 1;
            }

        };

        TNL::Containers::StaticArray<3, int> begin1{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end1{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceTypeHost>(begin1, end1, meshPoly);
         */

    }


    std::vector<boundaryConditionOutlet> boundary_vector_outlet_individual;
    std::vector<boundaryConditionInlet> boundary_vector_inlet_individual;
    std::vector<boundaryConditionInlet> boundary_vector_inlet;
    std::vector<boundaryConditionOutlet> boundary_vector_outlet;

    boundaryConditionInlet boundary_condition_inlet;

    LBMConstantsPointer Constants;
    LBMDataPointer Data;
};


#endif //GEOMETRYMESHER_H

#pragma clang diagnostic pop