#include "../traits/LBMTraits.h"

class LBMData
{
public:

    using DeviceType = LBMTraits::DeviceType;
    using DeviceTypeHost = LBMTraits::DeviceTypeHost;

    using ArrayTypeMesh = typename LBMTraits::ArrayTypeMesh;
    using ArrayTypeMeshHost = typename LBMTraits::ArrayTypeMeshHost;

    using ArrayTypeBoundaryInlet = typename LBMTraits::ArrayTypeBoundaryInlet;
    using ArrayTypeBoundaryInletHost = typename LBMTraits::ArrayTypeBoundaryInletHost;

    using ArrayTypeBoundaryOutlet = typename LBMTraits::ArrayTypeBoundaryOutlet;
    using ArrayTypeBoundaryOutletHost = typename LBMTraits::ArrayTypeBoundaryOutletHost;

    using ArrayTypeBoundaryWall = typename LBMTraits::ArrayTypeBoundaryWall;
    using ArrayTypeBoundaryWallHost = typename LBMTraits::ArrayTypeBoundaryWallHost;

    using ArrayTypeVariablesScalar = typename LBMTraits::ArrayTypeVariablesScalar;

    using ArrayTypeVariablesVector = typename LBMTraits::ArrayTypeVariablesVector;

    using ArrayTypeDFunction = typename LBMTraits::ArrayTypeDFunction;

    //mesh - host
    ArrayTypeMeshHost meshFluidHost;
    ArrayTypeBoundaryInletHost meshBoundaryInletHost;
    ArrayTypeBoundaryOutletHost meshBoundaryOutletHost;
    ArrayTypeBoundaryWallHost meshBoundaryWallHost;

    //mesh - device
    ArrayTypeMesh meshFluid;
    ArrayTypeBoundaryInlet meshBoundaryInlet;
    ArrayTypeBoundaryOutlet meshBoundaryOutlet;
    ArrayTypeBoundaryWall meshBoundaryWall;

    // variables
    ArrayTypeVariablesScalar rho;
    ArrayTypeVariablesVector u;
    ArrayTypeVariablesVector u_error;


    // distribution function

    ArrayTypeDFunction df;
    ArrayTypeDFunction df_post;

};