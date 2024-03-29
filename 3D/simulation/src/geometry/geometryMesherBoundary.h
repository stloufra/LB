#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

#ifndef GEOMETRYMESHER_H
#define GEOMETRYMESHER_H

#pragma once

#include <string>
#include <omp.h>


#include <vector>
#include "geometryObjectCuboid.h"
#include "../traits/LBMTraits.h"
#include "../data/LBMData.h"
#include "../data/LBMConstants.h"

#include <TNL/Algorithms/parallelFor.h>

using namespace TNL;
using namespace TNL::Algorithms;

class geometryMesherBoundary {
public:

    using RealType = typename LBMTraits::RealType;
    using DeviceType = typename LBMTraits::DeviceType;
    using DeviceTypeHost = typename LBMTraits::DeviceTypeHost;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;
    using ArrayTypeMeshHost = typename LBMTraits::ArrayTypeMeshHost;
    using ArrayTypeBoundaryInletHost = typename LBMTraits::ArrayTypeBoundaryInletHost;
    using ArrayTypeBoundaryOutletHost = typename LBMTraits::ArrayTypeBoundaryOutletHost;


    geometryMesherBoundary() = delete;

    geometryMesherBoundary(LBMConstantsPointer Constants_,
                           LBMDataPointer Data_) {

        Constants = Constants_;
        Data = Data_;
    }


    void meshingBoundaryConditionInlet(geometryObjectCuboid &cuboid, const Vector& normal_in, const Vector& velocity_in, bool verbose) {
        double lookx, looky, lookz;


        int num = 0;


        for (int k = 0; k < Constants->dimZ_int; k++) {

            for (int j = 0; j < Constants->dimY_int; j++) {

                for (int i = 0; i < Constants->dimX_int; i++) {

                    lookx = (static_cast<double>(i) / Constants->resolution_factor + Constants->BBminx -
                             Constants->additional_factor);
                    looky = (static_cast<double>(j) / Constants->resolution_factor + Constants->BBminy -
                             Constants->additional_factor);
                    lookz = (static_cast<double>(k) / Constants->resolution_factor + Constants->BBminz -
                             Constants->additional_factor);

                    d3 pnt_ask = {lookx, looky, lookz};

                    if (cuboid.isInside(pnt_ask) && Data->meshFluidHost(i, j, k) != 0) {
                        Data->meshFluidHost(i, j, k) = cuboid.id;


                        boundary_condition_inlet.normal = normal_in;
                        boundary_condition_inlet.velocity = velocity_in;
                        boundary_condition_inlet.vertex = {i, j, k};

                        boundary_vector_inlet.push_back(boundary_condition_inlet);
                        num += 1;
                    }
                }
            }
        }


        if (verbose) {
            std::cout << "\nMeshed cuboid id " << cuboid.id << " as Inlet.\n";
            std::cout << "Boundary vertexes number: " << num << std::endl;
        }
    }

    void compileBoundaryArrayInlets(bool verbose) {
        Constants->inlet_num = boundary_vector_inlet.size();
        std::cout << Constants->inlet_num << std::endl;

        Data->meshBoundaryInletHost.setSizes(Constants ->inlet_num);


        int i = 0;
        for (boundaryConditionInlet BC: boundary_vector_inlet) {
            Data->meshBoundaryInletHost[i] = BC;
            i += 1;
        }

        if (verbose) {
            std::cout << "\nCreated boundary Array Inlets " << ".\n";
            std::cout << "Boundary vertexes number: " << Constants -> inlet_num << std::endl;
        }

    }

    void meshingBoundaryConditionOutlet(geometryObjectCuboid &cuboid, const Vector& normal_in, RealType density, bool verbose) {
        double lookx, looky, lookz;

        boundaryConditionOutlet boundary_condition_outlet;

        int num = 0;


        for (int k = 0; k < Constants->dimZ_int; k++) {

            for (int j = 0; j < Constants->dimY_int; j++) {

                for (int i = 0; i < Constants->dimX_int; i++) {

                    lookx = (static_cast<double>(i) / Constants->resolution_factor + Constants->BBminx -
                             Constants->additional_factor);
                    looky = (static_cast<double>(j) / Constants->resolution_factor + Constants->BBminy -
                             Constants->additional_factor);
                    lookz = (static_cast<double>(k) / Constants->resolution_factor + Constants->BBminz -
                             Constants->additional_factor);

                    d3 pnt_ask = {lookx, looky, lookz};

                    if (cuboid.isInside(pnt_ask) && Data->meshFluidHost(i, j, k) != 0) {
                        Data->meshFluidHost(i, j, k) = cuboid.id;

                        boundary_condition_outlet.normal = normal_in;
                        boundary_condition_outlet.density = density;
                        boundary_condition_outlet.vertex = {i, j, k};

                        boundary_vector_outlet.push_back(boundary_condition_outlet);
                        num += 1;
                    }
                }
            }
        }

        /*for (boundaryConditionOutlet BC: boundary_vector_outlet_individual) {
            boundary_vector_outlet.push_back(BC);
        }*/

        if (verbose) {
            std::cout << "\nMeshed cuboid id " << cuboid.id << " as outlet.\n";
            std::cout << "Boundary vertexes number: " << num << std::endl;
        }
    }

    void compileBoundaryArrayOutlets(bool verbose) {
        Data->meshBoundaryOutletHost.setSizes(boundary_vector_outlet.size());
        Constants->outlet_num = boundary_vector_outlet.size();

        int j = 0;
        for (boundaryConditionOutlet BC: boundary_vector_outlet) {
            Data->meshBoundaryOutletHost[j] = BC;
            j += 1;
        }

        if (verbose) {
            std::cout << "\nCreated boundary Array Outlet.\n";
            std::cout << "Boundary vertexes number: " << boundary_vector_outlet.size() << std::endl;
        }
    }

    void meshingBoundaryWall(bool verbose) {

        boundaryConditionWall boundary_condition_wall;

        std::vector<boundaryConditionWall> boundary_vector_wall;

        for (int k = 0; k < Constants->dimZ_int; k++) {

            for (int j = 0; j < Constants->dimY_int; j++) {

                for (int i = 0; i < Constants->dimX_int; i++) {

                    if (Data->meshFluidHost(i, j, k) > 0) {
                        for (int l: {-1, 0, 1}) {
                            for (int m: {-1, 0, 1}) {
                                for (int n: {-1, 0, 1}) {
                                    if (i + l >= 0 && j + m >= 0 && k + n >= 0 && i + l < Constants->dimX_int &&
                                        j + m < Constants->dimY_int && k + n < Constants->dimZ_int) {
                                        if (Data->meshFluidHost(i + l, j + m, k + n) == 0) {
                                            Data->meshFluidHost(i, j, k) = 2;
                                            boundary_condition_wall.vertex = {i, j, k};
                                            boundary_vector_wall.push_back(boundary_condition_wall);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        Data->meshBoundaryWallHost.setSizes(boundary_vector_wall.size());
        Constants->wall_num = boundary_vector_wall.size();

        int k = 0;

        for (boundaryConditionWall BC: boundary_vector_wall) {
            Data->meshBoundaryWallHost[k] = BC;
            k += 1;
        }

        if (verbose) {
            std::cout << "\nMeshed and created boundary Array Wall.\n";
            std::cout << "Wall vertexes number: " << boundary_vector_wall.size() << std::endl;
        }


    }

    void arrayTransfer(bool verbose) {
        Data->meshFluid = Data->meshFluidHost;
        Data->meshBoundaryWall = Data->meshBoundaryWallHost;
        Data->meshBoundaryInlet = Data->meshBoundaryInletHost;
        Data->meshBoundaryOutlet = Data->meshBoundaryOutletHost;

        if (verbose) {
            if (std::is_same_v<DeviceType, TNL::Devices::Cuda>) {
                std::cout << "\n Trasfered Arrays From Host to Cuda.\n";
            } else if (std::is_same_v<DeviceType, TNL::Devices::Host>) {
                std::cout << "\n Trasfered Arrays From Host to Host.\n";
            }
        }
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