#include "../traits/LBMTraits.h"
#include <string>

struct LBMConstants {
    using RealType = typename LBMTraits::RealType;
    using VectorType = typename LBMTraits::VectorType;
    using VectorTypeInt = typename LBMTraits::VectorTypeInt;
    using string = typename std::string;

    //mesh
    RealType BBminx;                //[m]
    RealType BBmaxx;                //[m]
    RealType BBminy;                //[m]
    RealType BBmaxy;                //[m]
    RealType BBminz;                //[m]
    RealType BBmaxz;                //[m]

    int dimX_int;                   //[1]
    int dimY_int;                   //[1]
    int dimZ_int;                   //[1]

    RealType resolution_factor;     //[1]
    RealType additional_factor;     //[1]

    d3 point_outside;               //[m]

    string file_name;               //[word]

    int inlet_num;                  //[1]
    int outlet_num;                 //[1]
    int symmetry_num;               //[1]
    int periodic_num;               //[1]
    int periodicDP_num;             //[1]
    int wall_num;                   //[1]


    //physical
    RealType ny_fyz;
    RealType u_guess_fyz;           //[m/s]
    RealType L_fyz;                 //[m]
    RealType rho_fyz;               //[kg/m3]
    RealType Re;                    //[1]

    RealType Fx_fyz;                //[kg/m3/s2]
    RealType Fy_fyz;                //[kg/m3/s2]
    RealType Fz_fyz;                //[kg/m3/s2]
    RealType conversion_factor_fyz; //[1]  converting to meters 1/1000 mm; 1/100 cm; etc.

    //lattice
    RealType tau;
    RealType ny;
    RealType omega;
    RealType Fx;
    RealType Fy;
    RealType g;

    RealType U_lb;                  //[1]
    RealType U_inf = 0;                  //[1]
    RealType u_max;
    RealType Cl;
    RealType Ct;
    RealType Cm;
    RealType Cu_inverse;
    RealType Crho;
    RealType Cpressure;
    RealType Cu;

    //initialization
    VectorType VelocityInit = (0.f, 0.f, 0.f); //default
    string InitFileName;
    //turbulence
    const RealType CLES = 0.173f; //0.173f;

    //boundary dumping

    RealType omegaDumpingLow = -1.f; // if not set in main, dumping not performed

    RealType omegaDumpingHigh = -1.f; // if not set in main, dumping not performed
    //simulation
    RealType err;


    RealType time;                      //[s]
    RealType plot_every = -1.f;         //[s]
    RealType err_every = -1.f;          //[s]
    RealType MomentAvg_every = -1.f;     //[s]
    RealType MomentAvgStart = -1.f;     //[s]

    int plot_every_it = -1;
    int err_every_it = -1;

    int iterationsMomentAvg = -1;        // number of iteration to be time averaged
    int iterationsMomentAvgStart= -1;   // number of iteration to be time averaged

    int TimeAvgCounter = 0;         // number of iterations to be time averaged counter
    bool timeAveraged;              // if true, already time averaged

    int iterations;
    int probe_every_it = -1.f;             // number of iterations for probe sampling
    int probe_iterations = -1.f;
    VectorType ProbeLocation = (0.f, 0.f, 0.f); //default
    VectorTypeInt ProbeLocationLat = (0, 0, 0); //default



    int Nvel;
    const RealType cs = 0.57735f; //1.f/sqrt(3.f);
    const RealType cs2 = 1.f/3.f;
    const RealType cs4 = 1.f/9.f;

    enum PATCH {
        solid = 0,
        fluid = 1,
        wall = 2,
        //inlet outlet >2
    };
};