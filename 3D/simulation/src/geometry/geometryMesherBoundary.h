#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

#ifndef GEOMETRYMESHER_H
#define GEOMETRYMESHER_H

#pragma once

#include <string>
#include <omp.h>

#include <vector>
#include <tuple>
#include <map>
#include <unordered_map>
#include <cmath>
#include "geometryObjectCuboid.h"
#include "../traits/LBMTraits.h"
#include "../data/LBMData.h"
#include "../data/LBMConstants.h"

#include <TNL/Algorithms/parallelFor.h>

using namespace TNL;
using namespace TNL::Algorithms;

class geometryMesherBoundary
{
public:
    using RealType = typename LBMTraits::RealType;
    using VectorType = typename LBMTraits::VectorType;
    using DeviceType = typename LBMTraits::DeviceType;
    using DeviceTypeHost = typename LBMTraits::DeviceTypeHost;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;
    using ArrayTypeMeshHost = typename LBMTraits::ArrayTypeMeshHost;
    using ArrayTypeBoundaryInletHost = typename LBMTraits::ArrayTypeBoundaryInletHost;
    using ArrayTypeBoundaryOutletHost = typename LBMTraits::ArrayTypeBoundaryOutletHost;
    using PeriodicMap = typename LBMTraits::PeriodicMap;
    using PeriodicMapHost = typename LBMTraits::PeriodicMapHost;
    using PeriodicMapDevice = typename LBMTraits::PeriodicMapDevice;


    geometryMesherBoundary() = delete;

    geometryMesherBoundary(LBMConstantsPointer Constants_,
                           LBMDataPointer Data_)
    {
        Constants = Constants_;
        Data = Data_;
    }

    void meshingBoundaryConditionInletUniform(geometryObjectCuboid& cuboid,
                                              const Vector& normal_out,
                                              const Vector& velocity_in,
                                              bool verbose)
    {
        if (velocity_in.x() == velocity_in.y() && velocity_in.x() == velocity_in.z())
        {
            throw std::invalid_argument("Expected a 3D vector for velocity_in. If valid comment me.");
        }
        // Function to fill the inlet velocities with uniform velocity profile.
        // So far supports only inlets which are perpendicular to one of the axis.

        if (std::abs(normal_out.x() + normal_out.y() + normal_out.z()) != 1)
        {
            printf("Your normal is not tangent to one of the axes.\n");
        }

        RealType lookx, looky, lookz;


        int num = 0;


        for (int k = 0; k < Constants->dimZ_int; k++)
        {
            for (int j = 0; j < Constants->dimY_int; j++)
            {
                for (int i = 0; i < Constants->dimX_int; i++)
                {
                    lookx = (static_cast<double>(i) + 0.5 - Constants->additional_factor) / Constants->resolution_factor
                        + Constants->BBminx; // 0.5 for halfway between wall and first node
                    looky = (static_cast<double>(j) + 0.5 - Constants->additional_factor) / Constants->resolution_factor
                        + Constants->BBminy;
                    lookz = (static_cast<double>(k) + 0.5 - Constants->additional_factor) / Constants->resolution_factor
                        + Constants->BBminz;

                    d3 pnt_ask = {lookx, looky, lookz};

                    if (cuboid.isInside(pnt_ask) && Data->meshFluidHost(i, j, k) > 1)
                    {
                        if (neighboursFluidSearch(pnt_ask, 0))
                        {
                            boundary_condition_inlet.normal = normal_out;
                            boundary_condition_inlet.velocity = velocity_in;
                            boundary_condition_inlet.vertex = {i, j, k};

                            if (Data->meshFluidHost(i, j, k) > 10)
                            {
                                // count begins from +1
                                boundary_condition_inlet.regular = 0;
                            }
                            else
                            {
                                boundary_condition_inlet.regular = 1;
                            }

                            Data->meshFluidHost(i, j, k) = cuboid.id; // - Data->meshFluidHost(i, j, k);

                            boundary_vector_inlet.push_back(boundary_condition_inlet);
                            num += 1;
                        }
                    }
                }
            }
        }


        if (verbose)
        {
            std::cout << "\nMeshed cuboid id " << cuboid.id << " as Inlet.\n";
            std::cout << "Boundary vertexes number: " << num << std::endl;
        }
    }


    void meshingBoundaryConditionInletParaboloidCircle(geometryObjectCuboid& cuboid, // location of inlet
                                                       const d3& inletCenter, // center of inlet
                                                       const RealType& inletRadius, // inlet radius
                                                       const Vector& normal_out,
                                                       const RealType& mean_velocity_in,
                                                       bool verbose)
    {
        // Function to fill the inlet velocities with parabolic velocity profile.
        // So far supports only inlets which are perpendicular to one of the axis.

        if (std::abs(normal_out.x() + normal_out.y() + normal_out.z()) != 1)
        {
            printf("Your normal is not tangent to one of the axes.\n");
        }


        RealType lookx, looky, lookz;
        VectorType Velocity(0.f, 0.f, 0.f);


        meshingBoundaryConditionInletUniform(cuboid, normal_out, Velocity, 0);

        int num = 0;

        for (boundaryConditionInlet& BC : boundary_vector_inlet)
        {
            Vertex vert = BC.vertex;

            lookx = (static_cast<double>(vert.x) + 0.5 - Constants->additional_factor) / Constants->resolution_factor +
                Constants->BBminx; // 0.5 for halfway between wall and first node
            looky = (static_cast<double>(vert.y) + 0.5 - Constants->additional_factor) / Constants->resolution_factor +
                Constants->BBminy;
            lookz = (static_cast<double>(vert.z) + 0.5 - Constants->additional_factor) / Constants->resolution_factor +
                Constants->BBminz;

            if (cuboid.isInside({lookx, looky, lookz}))
            {
                // Calculate the distance from the inlet center
                RealType distance = std::sqrt(
                    std::pow(lookx - inletCenter.x, 2) +
                    std::pow(looky - inletCenter.y, 2) +
                    std::pow(lookz - inletCenter.z, 2)
                );

                // Check if the point is within the inlet radius
                if (distance <= inletRadius)
                {
                    // Calculate the parabolic velocity profile
                    RealType parabolicVelocity = 3 * mean_velocity_in * (1.0 - std::pow(distance / inletRadius, 2));

                    if (parabolicVelocity < 0)
                    {
                        parabolicVelocity = 0;
                    }


                    // Set the velocity components based on the orientation vector
                    BC.velocity.x() = -parabolicVelocity * normal_out.x();
                    BC.velocity.y() = -parabolicVelocity * normal_out.y();
                    BC.velocity.z() = -parabolicVelocity * normal_out.z();

                    num++;
                }
            }
        }


        if (verbose)
        {
            std::cout << "\nMeshed cuboid id " << cuboid.id << " as Inlet.\n";
            std::cout << "Boundary vertexes number: " << num << std::endl;
        }
    }

    void meshingBoundaryConditionInletParaboloidRectangle(geometryObjectCuboid& cuboid,
                                                          const d3& inletCenter,
                                                          RealType& inletDimX,
                                                          RealType& inletDimY,
                                                          RealType& inletDimZ,
                                                          const Vector& normal_out,
                                                          const RealType& mean_velocity_in,
                                                          bool verbose)
    {
        // Function to fill the inlet velocities with parabolic velocity profile.
        // So far supports only inlets which are perpendicular to one of the axes.

        if (std::abs(normal_out.x() + normal_out.y() + normal_out.z()) != 1)
        {
            printf("Your normal is not tangent to one of the axes.\n");
        }

        int num = 0;

        RealType lookx, looky, lookz;

        VectorType Velocity(mean_velocity_in, mean_velocity_in, mean_velocity_in);
        Velocity = -(Velocity * normal_out); //element wise

        meshingBoundaryConditionInletUniform(cuboid, normal_out, Velocity, 0);


        for (boundaryConditionInlet& BC : boundary_vector_inlet)
        {
            Vertex vert = BC.vertex;

            lookx = (static_cast<double>(vert.x) + 0.5 - Constants->additional_factor) / Constants->resolution_factor +
                Constants->BBminx; // 0.5 for halfway between wall and first node
            looky = (static_cast<double>(vert.y) + 0.5 - Constants->additional_factor) / Constants->resolution_factor +
                Constants->BBminy;
            lookz = (static_cast<double>(vert.z) + 0.5 - Constants->additional_factor) / Constants->resolution_factor +
                Constants->BBminz;

            if (cuboid.isInside({lookx, looky, lookz}))
            {
                // Calculate the distance from the inlet center

                VectorType RevNormal(1 - std::abs(normal_out.x()), 1 - std::abs(normal_out.y()),
                                     1 - std::abs(normal_out.z()));


                //distances from corner with min x, min y and min z values

                RealType distance_x = lookx - (inletCenter.x - inletDimX / 2);
                RealType distance_y = looky - (inletCenter.y - inletDimY / 2);
                RealType distance_z = lookz - (inletCenter.z - inletDimZ / 2);


                // Calculate the parabolic approximate velocity profile
                // u(x,y)=2 Uavg x (LX-x)*((y (LY-y))/((((LX)/(2)))^(2) (((LY)/(2)))^(2)))

                RealType Part_x, Part_y, Part_z;

                if (normal_out.x() == 0)
                {
                    Part_x = distance_x * (inletDimX - distance_x);
                }
                else
                {
                    Part_x = 1.f;
                    inletDimX = 1.f;
                }

                if (normal_out.y() == 0)
                {
                    Part_y = distance_y * (inletDimY - distance_y);
                }
                else
                {
                    Part_y = 1.f;
                    inletDimY = 1.f;
                }

                if (normal_out.z() == 0)
                {
                    Part_z = distance_z * (inletDimZ - distance_z);
                }
                else
                {
                    Part_z = 1.f;
                    inletDimZ = 1.f;
                }

                RealType inletConstant = inletDimX * inletDimY * inletDimZ * inletDimX * inletDimY * inletDimZ / 16;

                RealType parabolicVelocity = 2 * mean_velocity_in * Part_x * Part_y * Part_z /
                    inletConstant;


                if (parabolicVelocity < 0)
                {
                    parabolicVelocity = 0;
                }

                // Set the velocity components based on the orientation vector
                BC.velocity.x() = -parabolicVelocity * normal_out.x();
                BC.velocity.y() = -parabolicVelocity * normal_out.y();
                BC.velocity.z() = -parabolicVelocity * normal_out.z();
                num++;
            }
        }

        if (verbose)
        {
            std::cout << "\nMeshed cuboid id " << cuboid.id << " as Inlet.\n";
            std::cout << "Boundary vertexes number: " << num << std::endl;
        }
    }

    void meshingBoundaryConditionInletVariableExponent(geometryObjectCuboid& cuboid,
                                                       const d3& inletCenter,
                                                       RealType& inletDimX,
                                                       RealType& inletDimY,
                                                       RealType& inletDimZ,
                                                       const Vector& normal_out,
                                                       const RealType& mean_velocity_in,
                                                       const RealType& exponent_n,
                                                       char measurement_axis,
                                                       bool verbose)
    {
        if ((measurement_axis == 'x' && normal_out.x() != 0) ||
            (measurement_axis == 'y' && normal_out.y() != 0) ||
            (measurement_axis == 'z' && normal_out.z() != 0))
        {
            printf("Error: Measurement axis cannot be the same as the normal direction.\n");
            return;
        }

        int num = 0;
        RealType lookx, looky, lookz;

        VectorType Velocity(mean_velocity_in, mean_velocity_in, mean_velocity_in);
        Velocity = -(Velocity * normal_out); //element wise

        meshingBoundaryConditionInletUniform(cuboid, normal_out, Velocity, 0);

        for (boundaryConditionInlet& BC : boundary_vector_inlet)
        {
            Vertex vert = BC.vertex;

            lookx = (static_cast<double>(vert.x) + 0.5 - Constants->additional_factor) / Constants->resolution_factor +
                Constants->BBminx;
            looky = (static_cast<double>(vert.y) + 0.5 - Constants->additional_factor) / Constants->resolution_factor +
                Constants->BBminy;
            lookz = (static_cast<double>(vert.z) + 0.5 - Constants->additional_factor) / Constants->resolution_factor +
                Constants->BBminz;

            if (cuboid.isInside({lookx, looky, lookz}))
            {
                RealType r = 0.0f;
                RealType radius = 0.f;

                if (measurement_axis == 'x')
                {
                    r = std::abs(lookx - inletCenter.x);
                    radius = inletDimX / 2.f;
                }
                else if (measurement_axis == 'y')
                {
                    r = std::abs(looky - inletCenter.y);
                    radius = inletDimY / 2.f;
                }
                else if (measurement_axis == 'z')
                {
                    r = std::abs(lookz - inletCenter.z);
                    radius = inletDimZ / 2.f;
                }


                RealType Umax = mean_velocity_in * (exponent_n + 1.f) / (exponent_n);

                RealType velocity = Umax * std::pow((1.0f - r / radius), 1.0f / exponent_n);

                //printf("My vel %f", velocity);

                if (velocity < 0)
                {
                    velocity = 0;
                }


                BC.velocity.x() = -velocity * normal_out.x();
                BC.velocity.y() = -velocity * normal_out.y();
                BC.velocity.z() = -velocity * normal_out.z();
                num++;
            }
        }

        if (verbose)
        {
            std::cout << "\nMeshed cuboid id " << cuboid.id << " as Inlet with variable exponent.\n";
            std::cout << "Boundary vertexes number: " << num << std::endl;
        }
    }

    void meshingBoundaryConditionOutlet(geometryObjectCuboid& cuboid,
                                        const Vector& normal_out,
                                        RealType density,
                                        bool verbose)
    {
        double lookx, looky, lookz;

        boundaryConditionOutlet boundary_condition_outlet;

        int num = 0;


        for (int k = 0; k < Constants->dimZ_int; k++)
        {
            for (int j = 0; j < Constants->dimY_int; j++)
            {
                for (int i = 0; i < Constants->dimX_int; i++)
                {
                    lookx = (static_cast<double>(i) + 0.5 - Constants->additional_factor) / Constants->resolution_factor
                        + Constants->BBminx; // 0.5 for halfway between wall and first node
                    looky = (static_cast<double>(j) + 0.5 - Constants->additional_factor) / Constants->resolution_factor
                        + Constants->BBminy;
                    lookz = (static_cast<double>(k) + 0.5 - Constants->additional_factor) / Constants->resolution_factor
                        + Constants->BBminz;


                    d3 pnt_ask = {lookx, looky, lookz};

                    if (cuboid.isInside(pnt_ask) && Data->meshFluidHost(i, j, k) > 1)
                    {
                        boundary_condition_outlet.normal = normal_out;
                        boundary_condition_outlet.density = density;
                        boundary_condition_outlet.vertex = {i, j, k};

                        if (Data->meshFluidHost(i, j, k) == 10)
                        {
                            // count begins from +1
                            boundary_condition_outlet.regular = 1;
                        }
                        else
                        {
                            boundary_condition_outlet.regular = 0;
                        }

                        Data->meshFluidHost(i, j, k) = cuboid.id; // -Data->meshFluidHost(i, j, k)


                        boundary_vector_outlet.push_back(boundary_condition_outlet);
                        num += 1;
                    }
                }
            }
        }


        if (verbose)
        {
            std::cout << "\nMeshed cuboid id " << cuboid.id << " as outlet.\n";
            std::cout << "Boundary vertexes number: " << num << std::endl;
        }
    }

    void meshingBoundaryConditionSymmetry(  geometryObjectCuboid& cuboid,
                                            const Vector& normal_out,
                                            bool verbose)
    {
        // x = 0, y=1, z =2


        //#define SYMMETRY
        double lookx, looky, lookz;

        boundaryConditionSymmetry boundary_condition_symmetry;

        int num = 0;

        for (int k = 0; k < Constants->dimZ_int; k++)
        {
            for (int j = 0; j < Constants->dimY_int; j++)
            {
                for (int i = 0; i < Constants->dimX_int; i++)
                {
                    lookx = (static_cast<double>(i) + 0.5 - Constants->additional_factor) / Constants->resolution_factor
                        + Constants->BBminx; // 0.5 for halfway between wall and first node
                    looky = (static_cast<double>(j) + 0.5 - Constants->additional_factor) / Constants->resolution_factor
                        + Constants->BBminy;
                    lookz = (static_cast<double>(k) + 0.5 - Constants->additional_factor) / Constants->resolution_factor
                        + Constants->BBminz;


                    d3 pnt_ask = {lookx, looky, lookz};

                    if (cuboid.isInside(pnt_ask) && Data->meshFluidHost(i, j, k) != 0)
                    {
                        boundary_condition_symmetry.vertex = {i, j, k};
                        boundary_condition_symmetry.normal = normal_out;

                        if (Data->meshFluidHost(i, j, k) == 10)
                        {
                            // count begins from +1
                            boundary_condition_symmetry.regular = 1;
                        }
                        else
                        {
                            boundary_condition_symmetry.regular = 0;
                        }

                        Data->meshFluidHost(i, j, k) = cuboid.id; // Data->meshFluidHost(i, j, k)


                        boundary_vector_symmetry.push_back(boundary_condition_symmetry);
                        num += 1;
                    }
                }
            }
        }


        if (verbose)
        {
            std::cout << "\nMeshed cuboid id " << cuboid.id << " as Symmetry.\n";
            std::cout << "Boundary vertexes number: " << num << std::endl;
        }
    }

    void meshingBoundaryConditionPeriodic(  geometryObjectCuboid& cuboid,
                                            const Vector& normal_out,
                                            int periodicIndex_,
                                            bool verbose)
    {
        // x = 0, y=1, z =2

        double lookx, looky, lookz;

        boundaryConditionPeriodic boundary_condition_periodic;

        int num = 0;

        int ii, jj, kk;

        for (int k = 0; k < Constants->dimZ_int; k++)
        {
            for (int j = 0; j < Constants->dimY_int; j++)
            {
                for (int i = 0; i < Constants->dimX_int; i++)
                {
                    lookx = (static_cast<double>(i) + 0.5 - Constants->additional_factor) / Constants->resolution_factor
                        + Constants->BBminx; // 0.5 for halfway between wall and first node
                    looky = (static_cast<double>(j) + 0.5 - Constants->additional_factor) / Constants->resolution_factor
                        + Constants->BBminy;
                    lookz = (static_cast<double>(k) + 0.5 - Constants->additional_factor) / Constants->resolution_factor
                        + Constants->BBminz;


                    d3 pnt_ask = {lookx, looky, lookz};

                    if (cuboid.isInside(pnt_ask) && Data->meshFluidHost(i, j, k) != 0)
                    {
                        boundary_condition_periodic.vertex = {i, j, k};
                        boundary_condition_periodic.normal = normal_out;
                        boundary_condition_periodic.periodicIndex = periodicIndex_;;

                        if (Data->meshFluidHost(i, j, k) == 10)
                        {
                            // count begins from +1
                            boundary_condition_periodic.regular = 1;
                            //Data->meshFluidHost(i, j, k) = -3;
                        }
                        else
                        {
                            boundary_condition_periodic.regular = 0;
                            //Data->meshFluidHost(i, j, k) = -100;
                        }

                        Data->meshFluidHost(i, j, k) = cuboid.id; // Data->meshFluidHost(i, j, k)


                        boundary_vector_periodic.push_back(boundary_condition_periodic);
                        num += 1;

                        ii = i;
                        jj = j;
                        kk = k;
                    }
                }
            }
        }

        printf("\n My index is %d",
               abs(ii * (int)normal_out.x()) + abs(jj * (int)normal_out.y()) + abs(kk * (int)normal_out.z()));

        if (verbose)
        {
            std::cout << "\nMeshed cuboid id " << cuboid.id << " as Periodic.\n";
            std::cout << "Boundary vertexes number: " << num << std::endl;
        }
    }

    void meshingBoundaryConditionPeriodicDP(geometryObjectCuboid& cuboid,
                                            const Vector& normal_out,
                                            RealType& DeltaP,
                                            int periodicIndex_,
                                            bool verbose)
    {
        // x = 0, y=1, z =2

        double lookx, looky, lookz;

        boundaryConditionPeriodicDP boundary_condition_periodicDP;

        int num = 0;

        int ii, jj, kk;

        for (int k = 0; k < Constants->dimZ_int; k++)
        {
            for (int j = 0; j < Constants->dimY_int; j++)
            {
                for (int i = 0; i < Constants->dimX_int; i++)
                {
                    lookx = (static_cast<double>(i) + 0.5 - Constants->additional_factor) / Constants->resolution_factor
                        + Constants->BBminx; // 0.5 for halfway between wall and first node
                    looky = (static_cast<double>(j) + 0.5 - Constants->additional_factor) / Constants->resolution_factor
                        + Constants->BBminy;
                    lookz = (static_cast<double>(k) + 0.5 - Constants->additional_factor) / Constants->resolution_factor
                        + Constants->BBminz;


                    d3 pnt_ask = {lookx, looky, lookz};

                    if (cuboid.isInside(pnt_ask) && Data->meshFluidHost(i, j, k) != 0)
                    {
                        boundary_condition_periodicDP.vertex = {i, j, k};
                        boundary_condition_periodicDP.normal = normal_out;
                        boundary_condition_periodicDP.periodicIndex = periodicIndex_;
                        boundary_condition_periodicDP.DeltaRho = DeltaP;

                        if (Data->meshFluidHost(i, j, k) == 10)
                        {
                            // count begins from +1
                            boundary_condition_periodicDP.regular = 1;
                        }
                        else
                        {
                            boundary_condition_periodicDP.regular = 0;
                        }

                        Data->meshFluidHost(i, j, k) = cuboid.id; // Data->meshFluidHost(i, j, k)


                        boundary_vector_periodicDP.push_back(boundary_condition_periodicDP);
                        num += 1;

                        ii = i;
                        jj = j;
                        kk = k;
                    }
                }
            }
        }

        printf("\n My index is %d ",
               abs(ii * (int)normal_out.x()) + abs(jj * (int)normal_out.y()) + abs(kk * (int)normal_out.z()));
        printf("and my periodicNeighbour is %d", periodicIndex_);

        if (verbose)
        {
            std::cout << "\nMeshed cuboid id " << cuboid.id << " as Periodic with DP.\n";
            std::cout << "Boundary vertexes number: " << num << std::endl;
        }
    }

    void meshingBoundaryWall(bool verbose)
    {
        boundaryConditionWall boundary_condition_wall;


        for (int k = 0; k < Constants->dimZ_int; k++)
        {
            for (int j = 0; j < Constants->dimY_int; j++)
            {
                for (int i = 0; i < Constants->dimX_int; i++)
                {
                    if (Data->meshFluidHost(i, j, k) > 0)
                    {
                        for (int l : {-1, 0, 1})
                        {
                            for (int m : {-1, 0, 1})
                            {
                                for (int n : {-1, 0, 1})
                                {
                                    if (i + l >= 0 && j + m >= 0 && k + n >= 0 && i + l < Constants->dimX_int &&
                                        j + m < Constants->dimY_int && k + n < Constants->dimZ_int)
                                    {
                                        if (Data->meshFluidHost(i + l, j + m, k + n) == 0)
                                        {
                                            Data->meshFluidHost(i, j, k) += 1;
                                        }
                                    }
                                }
                            }
                        }

                        if (Data->meshFluidHost(i, j, k) > 2)
                        {
                            //TODO: before 2
                            boundary_condition_wall.vertex = {i, j, k};
                            boundary_condition_wall.regular = 1;

                            if (Data->meshFluidHost(i, j, k) == 10){
                                boundary_condition_wall.regular = 1;
                            }
                            else{
                                boundary_condition_wall.regular = 0;
                            }
                            boundary_vector_wall.push_back(boundary_condition_wall);
                        }
                    }
                }
            }
        }

#pragma omp parallel for
        for( boundaryConditionWall BC : boundary_vector_wall)
        {
            if( BC.regular == 1)
            {
                BC.normal = {0,0,0};

                auto x = BC.vertex.x;
                auto y = BC.vertex.y;
                auto z = BC.vertex.z;

                int check =0 ;
                for (int l : {-1, 1})
                {
                    if (Data->meshFluidHost(x + l, y,z) == 0)
                    {
                        BC.normal = {l,0,0};
                        check++;
                    }

                    if (Data->meshFluidHost(x,y+l,z) == 0)
                    {
                        BC.normal = {0,l,0};
                        check++;

                    }

                    if (Data->meshFluidHost(x, y,z+l) == 0)
                    {
                        BC.normal = {0,0,l};
                        check++;
//#define DEBUG
#ifdef DEBUG
                        Vector testNorm(0,0,-1);

                        if(BC.normal == testNorm)
                        {
                            Data->meshFluidHost(x,y,z) = 400;
                        }
#endif
                    }

                    BC.slipVelocity = {0.f,0.f,0.f};

                    if(check >1)
                    {
                        printf("FAIL more than one possible normal for wall node.\n");
                    }

                }

            }
        }

        if (verbose)
        {
            std::cout << "\nMeshed Wall.\n";
            std::cout << "Wall vertexes number: " << boundary_vector_wall.size() << std::endl;
        }
    }


    template <typename T>
    void fillSinglePeriodicMap(const T periodicVector, PeriodicMap& periodic_map)
    {
        for (const auto& boundary : periodicVector)
        {
            for (const auto& candidate : periodicVector)
            {
                if (boundary.vertex.x != candidate.vertex.x &&
                    boundary.vertex.y != candidate.vertex.y &&
                    boundary.vertex.z != candidate.vertex.z &&
                    boundary.normal == -candidate.normal)
                {
                    periodic_map[boundary.periodicIndex] = {candidate.periodicIndex, boundary.normal};
                    periodic_map[candidate.periodicIndex] = {boundary.periodicIndex, candidate.normal};
                    break;
                }
            }
        }
    }


    void tryin(bool verbose)
    {
        if (boundary_vector_periodic.size() > 0)
        {
            fillSinglePeriodicMap(boundary_vector_periodic, periodic_map);

            if (verbose)
            {
                std::cout << "\n";
                for (const auto& [key, value] : periodic_map)
                {
                    auto [partner_index, normal] = value;
                    std::cout << "Periodic Index " << key
                        << " -> Partner Index " << partner_index
                        << " with Normal (" << normal.x() << ", " << normal.y() << ", " << normal.z() << ")\n";
                }
            }
        }

        if (boundary_vector_periodicDP.size() > 0)
        {
            PeriodicMap periodic_mapDp;
            fillSinglePeriodicMap(boundary_vector_periodicDP, periodic_mapDp);
            if (verbose)
            {
                std::cout << "\n";
                for (const auto& [key, value] : periodic_mapDp)
                {
                    auto [partner_index, normal] = value;
                    std::cout << "PeriodicDP Index " << key
                        << " -> PartnerDP Index " << partner_index
                        << " with Normal (" << normal.x() << ", " << normal.y() << ", " << normal.z() << ")\n";
                }
            }
            periodic_map.insert(periodic_mapDp.begin(), periodic_mapDp.end());

        }

        PeriodicMapHost pmD(periodic_map.size());
        PeriodicMapDevice pmDD(periodic_map.size());

        int i=0;
        for (const auto& [key, value] : periodic_map)
        {
            auto [partner_index, normal] = value;
            pmD[i].key = key;
            pmD[i].partner_index = partner_index;
            pmD[i].partner_normal = normal;
            i++;
        }

        pmDD = pmD;
        auto pmDview = pmDD.getView();

        //printf("HEREEE %d",pmDview[0].key);

        auto init = [ = ] __cuda_callable__( int i ) mutable
        {
            //printf("HEREEE %d",i);
            //printf("Trying hash periodic boundary vector key - %d, value - %d \n", pmDview[i].key, pmDview[i].partner_index );
            printf("HEREEE %d",pmDview[0].key);
            printf("HEREEE %d",pmDview[1].key);
        };
        parallelFor< DeviceType >( 0, periodic_map.size(), init );


    }

    void compileBoundaryArrayInlets(bool verbose)
    {
        Constants->inlet_num = boundary_vector_inlet.size();

        Data->meshBoundaryInletHost.setSizes(Constants->inlet_num);

        int i = 0;
        for (boundaryConditionInlet BC : boundary_vector_inlet)
        {
            Data->meshBoundaryInletHost[i] = BC;
            i += 1;
        }

        if (verbose)
        {
            std::cout << "\nCreated boundary Array Inlets " << ".\n";
            std::cout << "Boundary vertexes number: " << Constants->inlet_num << std::endl;
        }
    }

    void compileBoundaryArrayOutlets(bool verbose)
    {
        Data->meshBoundaryOutletHost.setSizes(boundary_vector_outlet.size());
        Constants->outlet_num = boundary_vector_outlet.size();

        int j = 0;
        for (boundaryConditionOutlet BC : boundary_vector_outlet)
        {
            Data->meshBoundaryOutletHost[j] = BC;
            j += 1;
        }

        if (verbose)
        {
            std::cout << "\nCreated boundary Array Outlet.\n";
            std::cout << "Boundary vertexes number: " << boundary_vector_outlet.size() << std::endl;
        }
    }

    void compileBoundaryArrayWall(bool verbose)
    {
        auto it = boundary_vector_wall.begin();

        int er = 0;


#if __cplusplus >= 202002L
        er = std::erase_if(boundary_vector_wall, [&](const boundaryConditionWall& BC) {
        return Data->meshFluidHost(BC.vertex.x, BC.vertex.y, BC.vertex.z) < 0;
        });
#elif __cplusplus <= 202002L
        while (it != boundary_vector_wall.end())
        {
            const auto& Ver = it->vertex;

            if (Data->meshFluidHost(Ver.x, Ver.y, Ver.z) < 0)
            {
                it = boundary_vector_wall.erase(it);
                ++er;
            }
            else
            {
                ++it;
            }
        }
#endif
        Constants->wall_num = boundary_vector_wall.size();

        Data->meshBoundaryWallHost.setSizes(boundary_vector_wall.size());

        int k = 0;

        for (boundaryConditionWall BC : boundary_vector_wall)
        {
            Data->meshBoundaryWallHost[k] = BC;
            ++k;
        }

        if (verbose)
        {
            std::cout << "\nCreated boundary Array Wall.\n";
            std::cout << "Erased " << er << " elements.\n";
            std::cout << "Wall vertexes number: " << boundary_vector_wall.size() << std::endl;
        }
    };

    void compileBoundaryArraySymmetry(bool verbose)
    {
        Constants->symmetry_num = boundary_vector_symmetry.size();

        Data->meshBoundarySymmetryHost.setSizes(Constants->symmetry_num);

        int k = 0;

        for (const auto& BC : boundary_vector_symmetry)
        {
            Data->meshBoundarySymmetryHost[k] = BC;
            ++k;
        }

        if (verbose)
        {
            std::cout << "\nCreated boundary Array Symmetry.\n";
            std::cout << "Symmetry vertexes number: " << boundary_vector_wall.size() << std::endl;
        }
    };

    void compileBoundaryArrayPeriodic(bool verbose)
    {
        Constants->periodic_num = boundary_vector_periodic.size();

        Data->meshBoundaryPeriodicHost.setSizes(Constants->periodic_num);

        int k = 0;

        for (const auto& BC : boundary_vector_periodic)
        {
            Data->meshBoundaryPeriodicHost[k] = BC;
            ++k;
        }

        if (verbose)
        {
            std::cout << "\nCreated boundary Array Periodic.\n";
            std::cout << "Periodic vertexes number: " << boundary_vector_periodic.size() << std::endl;
        }
    };

    void compileBoundaryArrayPeriodicDP(bool verbose)
    {
        Constants->periodicDP_num = boundary_vector_periodicDP.size();

        Data->meshBoundaryPeriodicDPHost.setSizes(Constants->periodicDP_num);

        int k = 0;

        for (const auto& BC : boundary_vector_periodicDP)
        {
            Data->meshBoundaryPeriodicDPHost[k] = BC;
            ++k;
        }

        if (verbose)
        {
            std::cout << "\nCreated boundary Array Periodic with Delta p.\n";
            std::cout << "Periodic Dp vertexes number: " << boundary_vector_periodicDP.size() << std::endl;
        }
    };

    void arrayTransfer(bool verbose)
    {
        Data->meshFluid = Data->meshFluidHost;

        if (boundary_vector_wall.size())
        {
            compileBoundaryArrayWall(verbose);
            Data->meshBoundaryWall = Data->meshBoundaryWallHost;
        }

        if (boundary_vector_inlet.size())
        {
            compileBoundaryArrayInlets(verbose);
            Data->meshBoundaryInlet = Data->meshBoundaryInletHost;
        }

        if (boundary_vector_outlet.size())
        {
            compileBoundaryArrayOutlets(verbose);
            Data->meshBoundaryOutlet = Data->meshBoundaryOutletHost;
        }

        if (boundary_vector_symmetry.size())
        {
            compileBoundaryArraySymmetry(verbose);
            Data->meshBoundarySymmetry = Data->meshBoundarySymmetryHost;
        }


        if (boundary_vector_periodic.size())
        {
            compileBoundaryArrayPeriodic(verbose);
            Data->meshBoundaryPeriodic = Data->meshBoundaryPeriodicHost;
        }

        if (boundary_vector_periodicDP.size())
        {
            compileBoundaryArrayPeriodicDP(verbose);
            Data->meshBoundaryPeriodicDP = Data->meshBoundaryPeriodicDPHost;
        }

        if (verbose)
        {
            if (std::is_same_v<DeviceType, TNL::Devices::Cuda>)
            {
                std::cout << "\n Trasfered Arrays From Host to Cuda.\n";
            }
            else if (std::is_same_v<DeviceType, TNL::Devices::Host>)
            {
                std::cout << "\n Trasfered Arrays From Host to Host.\n";
            }
        }
    }

    // HELPER FUNCTIONS

    bool neighboursFluidSearch(d3 Point, int TypeSearched)
    {
        int nx, ny, nz;

        // Iterate through neighbors
        for (int di = -1; di <= 1; ++di)
        {
            for (int dj = -1; dj <= 1; ++dj)
            {
                for (int dk = -1; dk <= 1; ++dk)
                {
                    if (di == 0 && dj == 0 && dk == 0) // Skip the central point
                    {
                        continue;
                    }

                    // Calculate neighbor indices
                    nx = Point.x + di;
                    ny = Point.y + dj;
                    nz = Point.z + dk;

                    // Check if the neighbor is within array bounds
                    if (nx >= 0 && nx < Constants->dimX_int &&
                        ny >= 0 && ny < Constants->dimY_int &&
                        nz >= 0 && nz < Constants->dimZ_int)
                    {
                        // Check if the neighbor has a value of TypeSearched
                        if (Data->meshFluidHost(nx, ny, nz) == TypeSearched)
                        {
                            return true; // Found a neighbor with a value of TypeSearched
                        }
                    }
                }
            }
        }

        return true; // No neighbor with a value of TypeSearched found
    }

    std::vector<boundaryConditionInlet> boundary_vector_inlet;
    std::vector<boundaryConditionOutlet> boundary_vector_outlet;
    std::vector<boundaryConditionSymmetry> boundary_vector_symmetry;
    std::vector<boundaryConditionPeriodic> boundary_vector_periodic;
    std::vector<boundaryConditionPeriodicDP> boundary_vector_periodicDP;
    std::vector<boundaryConditionWall> boundary_vector_wall;

    boundaryConditionInlet boundary_condition_inlet;

    PeriodicMap periodic_map;

    LBMConstantsPointer Constants;
    LBMDataPointer Data;
};


#endif //GEOMETRYMESHER_H

#pragma clang diagnostic pop
