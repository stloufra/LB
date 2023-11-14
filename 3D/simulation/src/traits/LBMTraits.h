#ifndef LBMTraits_H
#define LBMTraits_H

#include <TNL/Containers/NDArray.h>
#include <TNL/Containers/StaticArray.h>
#include <TNL/Pointers/SharedPointer.h>
#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Algorithms/reduce.h>
#include <TNL/Containers/StaticVector.h>
#include <vector>

#pragma once
using namespace TNL;
using namespace TNL::Algorithms;

typedef struct {
    int x, y, z;
} Vertex;

typedef struct {
    int x, y, z;
} Normal;

//typedef struct {
//    double ux, uy, uz;
//} Velocity;

using Vector = TNL::Containers::StaticVector< 3, float >;

typedef struct {
    Vertex vertex;
    Vector normal;
    Vector velocity;
} boundaryConditionInlet;

typedef struct {
    Vertex vertex;
    Vector normal;
    float density;
} boundaryConditionOutlet;

typedef struct {
    Vertex vertex;
} boundaryConditionWall;

typedef struct {
    double x, y, z;
} d3;

class LBMTraits {

public:


    using RealType = float;

    using DeviceType = TNL::Devices::Cuda;

    using DeviceTypeHost = TNL::Devices::Host;
    using VectorType = TNL::Containers::StaticVector< 3, RealType >;

    // ---------------- MESH -----------------------
    using ArrayTypeMeshHost = TNL::Containers::NDArray<         int,
                                                                TNL::Containers::SizesHolder<int, 0, 0, 0>,
                                                                std::index_sequence<0, 1, 2>,
                                                                DeviceTypeHost>;

    using ArrayTypeMesh = TNL::Containers::NDArray<             int,
                                                                TNL::Containers::SizesHolder<int, 0, 0, 0>,
                                                                std::index_sequence<0, 1, 2>,
                                                                DeviceType>;

    // ---------------- INLET/OUTLET -------------------
    using ArrayTypeBoundaryInletHost = TNL::Containers::NDArray<boundaryConditionInlet,
                                                                TNL::Containers::SizesHolder<int, 0>,
                                                                std::index_sequence<0>,
                                                                DeviceTypeHost>;

    using ArrayTypeBoundaryOutletHost = TNL::Containers::NDArray<boundaryConditionOutlet,
                                                                TNL::Containers::SizesHolder<int, 0>,
                                                                std::index_sequence<0>,
                                                                DeviceTypeHost>;


    using ArrayTypeBoundaryInlet = TNL::Containers::NDArray <   boundaryConditionInlet,
                                                                TNL::Containers::SizesHolder<int, 0>,
                                                                std::index_sequence<0>,
                                                                DeviceType>;

    using ArrayTypeBoundaryOutlet = TNL::Containers::NDArray <  boundaryConditionOutlet,
                                                                TNL::Containers::SizesHolder<int, 0>,
                                                                std::index_sequence<0>,
                                                                DeviceType>;

    // -------------------- WALL --------------------------
    using ArrayTypeBoundaryWall = TNL::Containers::NDArray <    boundaryConditionWall ,
                                                                TNL::Containers::SizesHolder<int, 0>,
                                                                std::index_sequence<0>,
                                                                DeviceType>;

    using ArrayTypeBoundaryWallHost = TNL::Containers::NDArray< boundaryConditionWall ,
                                                                TNL::Containers::SizesHolder<int, 0>,
                                                                std::index_sequence<0>,
                                                                DeviceTypeHost>;
    // -------------------- VARIABLES --------------------------
    using ArrayTypeVariablesScalar = TNL::Containers::NDArray<  RealType,
                                                                TNL::Containers::SizesHolder<int, 0, 0, 0>,
                                                                std::index_sequence<0, 1, 2>,
                                                                DeviceType>;

    using ArrayTypeVariablesVector = TNL::Containers::NDArray<  VectorType,
                                                                TNL::Containers::SizesHolder<int, 0, 0, 0>,
                                                                std::index_sequence<0, 1, 2>,
                                                                DeviceType>;

    using ArrayTypeVariablesScalarHost = TNL::Containers::NDArray<  RealType,
                                                                TNL::Containers::SizesHolder<int, 0, 0, 0>,
                                                                std::index_sequence<0, 1, 2>,
                                                                DeviceTypeHost>;

    using ArrayTypeVariablesVectorHost = TNL::Containers::NDArray<  VectorType,
                                                                TNL::Containers::SizesHolder<int, 0, 0, 0>,
                                                                std::index_sequence<0, 1, 2>,
                                                                DeviceTypeHost>;

    // -------------------- DFUNCTION --------------------------
    using ArrayTypeDFunction = TNL::Containers::NDArray<        RealType,
                                                                TNL::Containers::SizesHolder<int, 0, 0, 0, 0>,
                                                                std::index_sequence<0, 1, 2, 3>,
                                                                DeviceType>;

    using ArrayTypeDFunctionHost = TNL::Containers::NDArray<    RealType,
                                                                TNL::Containers::SizesHolder<int, 0, 0, 0, 0>,
                                                                std::index_sequence<0, 1, 2, 3>,
                                                                DeviceTypeHost>;


};

#endif
