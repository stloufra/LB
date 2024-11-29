#ifndef LBMTraits_H
#define LBMTraits_H

#include <TNL/Containers/NDArray.h>
#include <TNL/Containers/StaticArray.h>
#include <TNL/Pointers/SharedPointer.h>
#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Algorithms/reduce.h>
#include <TNL/Containers/StaticVector.h>
#include <TNL/Containers/Vector.h>
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
    Vector normal;
    bool regular;
} boundaryConditionSymmetry;

typedef struct {
    Vertex vertex;
    Vector normal;
    int periodicIndex;
    bool regular;
} boundaryConditionPeriodic;

typedef struct {
    Vertex vertex;
    Vector normal;
    float DeltaRho;
    int periodicIndex;
    bool regular;
} boundaryConditionPeriodicDP;

typedef struct {
    Vertex vertex;
    bool regular;
} boundaryConditionWall;

typedef struct {
    int key;
    int partner_index;
    Vector partner_normal;
} periodicHash;

typedef struct {
    double x, y, z;
} d3;

class LBMTraits {

public:

    // Helper to detect if a type exists
    template <typename, typename = std::void_t<>>
    struct is_defined : std::false_type {};

    // Specialization if the type exists
    template <typename T>
    struct is_defined<T, std::void_t<decltype(std::declval<T>())>> : std::true_type {};


    using RealType = float;

    using DeviceType = TNL::Devices::Cuda;
    using DeviceTypeHost = TNL::Devices::Host;

    using VectorType = Vector;
    using VectorTypeInt = TNL::Containers::StaticVector< 3, int >;
    using VectorTypeProbeCuda = TNL::Containers::Vector< RealType, DeviceType >;
    using VectorTypeProbeHost = TNL::Containers::Vector< RealType, DeviceTypeHost >;

    using TensorType = TNL::Containers::StaticVector<3, TNL::Containers::StaticVector< 3, RealType >>;

    using PeriodicMap= std::unordered_map<int, std::pair<int, VectorType>>;
    using PeriodicMapHost= TNL::Containers::Array<periodicHash, DeviceTypeHost>;
    using PeriodicMapDevice = TNL::Containers::Array<periodicHash, DeviceType>;

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

    // -------------------- PERIODIC -------------------------

    using ArrayTypeBoundaryPeriodic = TNL::Containers::NDArray <    boundaryConditionPeriodic ,
                                                                TNL::Containers::SizesHolder<int, 0>,
                                                                std::index_sequence<0>,
                                                                DeviceType>;

    using ArrayTypeBoundaryPeriodicHost = TNL::Containers::NDArray< boundaryConditionPeriodic ,
                                                                TNL::Containers::SizesHolder<int, 0>,
                                                                std::index_sequence<0>,
                                                                DeviceTypeHost>;

    using ArrayTypeBoundaryPeriodicDP = TNL::Containers::NDArray <    boundaryConditionPeriodicDP ,
                                                               TNL::Containers::SizesHolder<int, 0>,
                                                               std::index_sequence<0>,
                                                               DeviceType>;

    using ArrayTypeBoundaryPeriodicDPHost = TNL::Containers::NDArray< boundaryConditionPeriodicDP ,
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
