#include "../traits/LBMTraits.h"
#include <string>

struct LBMConstants {
    using RealType = typename LBMTraits::RealType;
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
    string mesh_name;               //[word]
    int inlet_num;                  //[1]
    int outlet_num;                    //[1]
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
    RealType u_max;
    RealType Cl;
    RealType Ct;
    RealType Cm;
    RealType Cu_inverse;
    RealType Crho;
    RealType Cu;

    //simulation
    RealType err;
    RealType time;                  //[s]  
    RealType plot_every;            //[s]
    int plot_every_it;
    int iterations;

    //TODO: to model
    int Nvel;
    const RealType cs = 1/sqrt(3.f);
    const RealType cs2 = 1/3.f;


    enum PATCH {
        solid = 0,
        fluid = 1,
        wall = 2,
        //inlet outlet >2
    };
};