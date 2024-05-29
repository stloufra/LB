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
    bool regular;
} boundaryConditionInlet;

typedef struct {
    Vertex vertex;
    Vector normal;
    float density;
    bool regular;
} boundaryConditionOutlet;

typedef struct {
    Vertex vertex;
    int normal;
    bool regular;
} boundaryConditionSymmetry;

typedef struct {
    Vertex vertex;
    bool regular;
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
    using TensorType = TNL::Containers::StaticVector<3, TNL::Containers::StaticVector< 3, RealType >>;

    using NDArray3DSequenceType = std::index_sequence<2, 1, 0>;

    using NDArray4DSequenceType = std::index_sequence<3, 2, 1, 0>;

    // ---------------- MESH -----------------------
    using ArrayTypeMeshHost = TNL::Containers::NDArray<         int,
                                                                TNL::Containers::SizesHolder<int, 0, 0, 0>,
                                                                NDArray3DSequenceType,
                                                                DeviceTypeHost>;

    using ArrayTypeMesh = TNL::Containers::NDArray<             int,
                                                                TNL::Containers::SizesHolder<int, 0, 0, 0>,
                                                                NDArray3DSequenceType,
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
    // -------------------- SYMMETRY -------------------------

    using ArrayTypeBoundarySymmetry = TNL::Containers::NDArray <    boundaryConditionSymmetry ,
                                                                TNL::Containers::SizesHolder<int, 0>,
                                                                std::index_sequence<0>,
                                                                DeviceType>;

    using ArrayTypeBoundarySymmetryHost = TNL::Containers::NDArray< boundaryConditionSymmetry ,
                                                                TNL::Containers::SizesHolder<int, 0>,
                                                                std::index_sequence<0>,
                                                                DeviceTypeHost>;


    // -------------------- VARIABLES ------------------------
    using ArrayTypeVariablesScalar = TNL::Containers::NDArray<  RealType,
                                                                TNL::Containers::SizesHolder<int, 0, 0, 0>,
                                                                NDArray3DSequenceType,
                                                                DeviceType>;

    using ArrayTypeVariablesVector = TNL::Containers::NDArray<  VectorType,
                                                                TNL::Containers::SizesHolder<int, 0, 0, 0>,
                                                                NDArray3DSequenceType,
                                                                DeviceType>;

    using ArrayTypeVariablesScalarHost = TNL::Containers::NDArray<  RealType,
                                                                TNL::Containers::SizesHolder<int, 0, 0, 0>,
                                                                NDArray3DSequenceType,
                                                                DeviceTypeHost>;

    using ArrayTypeVariablesVectorHost = TNL::Containers::NDArray<  VectorType,
                                                                TNL::Containers::SizesHolder<int, 0, 0, 0>,
                                                                NDArray3DSequenceType,
                                                                DeviceTypeHost>;

    // -------------------- DFUNCTION --------------------------
    using ArrayTypeDFunction = TNL::Containers::NDArray<        RealType,
                                                                TNL::Containers::SizesHolder<int, 0, 0, 0, 0>,
                                                                NDArray4DSequenceType,
                                                                DeviceType>;

    using ArrayTypeDFunctionHost = TNL::Containers::NDArray<    RealType,
                                                                TNL::Containers::SizesHolder<int, 0, 0, 0, 0>,
                                                                NDArray4DSequenceType,
                                                                DeviceTypeHost>;


};

#endif
