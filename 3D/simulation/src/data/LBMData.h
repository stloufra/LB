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

    using ArrayTypeBoundarySymmetry = typename LBMTraits::ArrayTypeBoundarySymmetry;
    using ArrayTypeBoundarySymmetryHost = typename LBMTraits::ArrayTypeBoundarySymmetryHost;

    using ArrayTypeBoundaryPeriodic = typename LBMTraits::ArrayTypeBoundaryPeriodic;
    using ArrayTypeBoundaryPeriodicHost = typename LBMTraits::ArrayTypeBoundaryPeriodicHost;

    using ArrayTypeVariablesScalar = typename LBMTraits::ArrayTypeVariablesScalar;
    using ArrayTypeVariablesScalarHost = typename LBMTraits::ArrayTypeVariablesScalarHost;

    using ArrayTypeVariablesVector = typename LBMTraits::ArrayTypeVariablesVector;
    using ArrayTypeVariablesVectorHost = typename LBMTraits::ArrayTypeVariablesVectorHost;

    using ArrayTypeDFunction = typename LBMTraits::ArrayTypeDFunction;
    using ArrayTypeDFunctionHost = typename LBMTraits::ArrayTypeDFunctionHost;

    //mesh - host
    ArrayTypeMeshHost meshFluidHost;
    ArrayTypeBoundaryInletHost meshBoundaryInletHost;
    ArrayTypeBoundaryOutletHost meshBoundaryOutletHost;
    ArrayTypeBoundarySymmetryHost meshBoundarySymmetryHost;
    ArrayTypeBoundaryPeriodicHost meshBoundaryPeriodicHost;
    ArrayTypeBoundaryWallHost meshBoundaryWallHost;

    //mesh - device
    ArrayTypeMesh meshFluid;
    ArrayTypeBoundaryInlet meshBoundaryInlet;
    ArrayTypeBoundaryOutlet meshBoundaryOutlet;
    ArrayTypeBoundarySymmetry meshBoundarySymmetry;
    ArrayTypeBoundaryPeriodic meshBoundaryPeriodic;
    ArrayTypeBoundaryWall meshBoundaryWall;

    // variables
    ArrayTypeVariablesScalar rho;
    ArrayTypeVariablesScalar p;
    ArrayTypeVariablesVector u;

    ArrayTypeVariablesScalar rhoTimeAvg;
    ArrayTypeVariablesScalar pTimeAvg;
    ArrayTypeVariablesVector uTimeAvg;

    ArrayTypeVariablesVector u_error;
    ArrayTypeVariablesScalar omega;

    ArrayTypeVariablesScalarHost rho_out;
    ArrayTypeVariablesScalarHost omega_out;
    ArrayTypeVariablesScalarHost p_out;
    ArrayTypeVariablesVectorHost u_out;



    // distribution function

    ArrayTypeDFunction df;
    ArrayTypeDFunction df_post;

};