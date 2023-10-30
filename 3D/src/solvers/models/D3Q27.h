#ifndef D3Q27_H
#define  D3Q27_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>
#include <omp.h>
#include "../../geometry/geometryMesher.h"
#include "../../traits/LBMTraits.h"
#include <TNL/Containers/StaticArray.h>
#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Algorithms/reduce.h>

using namespace TNL;
using namespace TNL::Algorithms;

//template< typename INIT >//typename COL //TODO implement
class D3Q27 {
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;


public:

    static void initializeSim(LBMDataPointer &Data, LBMConstantsPointer &Constants, const Vector& u_0) {


        auto rho_view = Data->rho.getView();
        auto u_view = Data->u.getView();
        auto mesh_view = Data->meshFluid.getView();
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();

        const RealType Cl = Constants->Cl;
        const RealType Ct = Constants->Ct;
        const RealType Cm = Constants->Cm;
        const RealType Cu_inverse = Constants->Cu_inverse;
        const RealType Cu = Constants->Cu;
        const int Nvel = Constants->Nvel;


        //const int (&c)[27][3] = this->c; //cx_post previously
        //const RealType (&weight)[27] = weight;

        constexpr int c[27][3] = {
                {0,  0,  0},

                {1,  0,  0},
                {-1, 0,  0},
                {0,  0,  -1},
                {0,  0,  1},
                {0,  -1, 0},
                {0,  1,  0},

                {1,  0,  -1},
                {-1, 0,  1},
                {1,  0,  1},
                {-1, 0,  -1},
                {-1, -1, 0},
                {1,  1,  0},
                {0,  1,  -1},
                {0,  -1, 1},
                {-1, 1,  0},
                {1,  -1, 0},
                {0,  1,  1},
                {0,  -1, -1},

                {-1, 1,  -1},
                {1,  -1, 1},
                {-1, -1, 1},
                {1,  1,  -1},
                {1,  -1, -1},
                {-1, 1,  1},
                {-1, -1, -1},
                {1,  1,  1}};

        constexpr static RealType weight[27] = {8.f / 27, // 0
                                                2.f / 27, 2.f / 27, 2.f / 27, 2.f / 27, 2.f / 27, 2.f / 27, // 1-6
                                                1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54,
                                                1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, // 7- 18
                                                1.f / 216, 1.f / 216, 1.f / 216, 1.f / 216, 1.f / 216,
                                                1.f / 216, 1.f / 216, 1.f / 216}; // 19-26


        RealType rho_0 = Constants->rho_fyz;
        //VectorType u_0(0.25f, 0, 0); // TODO after solving fetch


        auto f_equilibrium = [=]
        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k,
        const int &vel ) mutable
        {
            RealType uc, u2;

            uc = c[vel][0] * u_view(i, j, k).x()
                 + c[vel][1] * u_view(i, j, k).y()
                 + c[vel][2] * u_view(i, j, k).z();

            u2 = u_view(i, j, k).x() * u_view(i, j, k).x()
                 + u_view(i, j, k).y() * u_view(i, j, k).y()
                 + u_view(i, j, k).z() * u_view(i, j, k).z();


            return weight[vel] * rho_view(i, j, k) * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };


        auto init_variables = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i ) mutable
        {
            if (mesh_view(i.x(), i.y(), i.z()) != 0) {
                rho_view(i.x(), i.y(), i.z()) = rho_0 / Cm * Cl * Cl * Cl;
                u_view(i.x(), i.y(), i.z()).x() = u_0.x() * Cu_inverse;
                u_view(i.x(), i.y(), i.z()).y() = u_0.y() * Cu_inverse;
                u_view(i.x(), i.y(), i.z()).z() = u_0.z() * Cu_inverse;

            } else {
                rho_view(i.x(), i.y(), i.z()) = 0.f;
                u_view(i.x(), i.y(), i.z()).x() = 0.f;
                u_view(i.x(), i.y(), i.z()).y() = 0.f;
                u_view(i.x(), i.y(), i.z()).z() = 0.f;
            }
        };

        auto init_df = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i ) mutable
        {
            for (int vel = 0; vel < Constants->Nvel; vel++) {
                if (mesh_view(i.x(), i.y(), i.z()) !=0) {
                    df_view(i.x(), i.y(), i.z(), vel) = f_equilibrium(i.x(), i.y(), i.z(), vel);
                    df_post_view(i.x(), i.y(), i.z(), vel) = f_equilibrium(i.x(), i.y(), i.z(), vel);

                } else {
                    df_view(i.x(), i.y(), i.z(), vel) = 0.f;
                    df_post_view(i.x(), i.y(), i.z(), vel) = 0.f;
                }
            }
        };


        TNL::Containers::StaticArray<3, int> begin1{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end1{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        // todo try mesh_view.getSizes()
        parallelFor<DeviceType>(begin1, end1, init_variables);

        parallelFor<DeviceType>(begin1, end1, init_df);
        //printf( "init_df \n");



    };

    static void collision(LBMDataPointer &Data, LBMConstantsPointer &Constants) {
        auto rho_view = Data->rho.getView();
        auto u_view = Data->u.getView();
        auto mesh_view = Data->meshFluid.getView();
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();



        constexpr int c[27][3] = {
                {0,  0,  0},

                {1,  0,  0},
                {-1, 0,  0},
                {0,  0,  -1},
                {0,  0,  1},
                {0,  -1, 0},
                {0,  1,  0},

                {1,  0,  -1},
                {-1, 0,  1},
                {1,  0,  1},
                {-1, 0,  -1},
                {-1, -1, 0},
                {1,  1,  0},
                {0,  1,  -1},
                {0,  -1, 1},
                {-1, 1,  0},
                {1,  -1, 0},
                {0,  1,  1},
                {0,  -1, -1},

                {-1, 1,  -1},
                {1,  -1, 1},
                {-1, -1, 1},
                {1,  1,  -1},
                {1,  -1, -1},
                {-1, 1,  1},
                {-1, -1, -1},
                {1,  1,  1}};

        constexpr static RealType weight[27] = {8.f / 27, // 0
                                                2.f / 27, 2.f / 27, 2.f / 27, 2.f / 27, 2.f / 27, 2.f / 27, // 1-6
                                                1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54,
                                                1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, // 7- 18
                                                1.f / 216, 1.f / 216, 1.f / 216, 1.f / 216, 1.f / 216,
                                                1.f / 216, 1.f / 216, 1.f / 216}; // 19-26

        auto f_equilibrium = [=]
        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k,
        const int &vel ) mutable
        {
            RealType uc, u2;

            uc = c[vel][0] * u_view(i, j, k).x()
                 + c[vel][1] * u_view(i, j, k).y()
                 + c[vel][2] * u_view(i, j, k).z();

            u2 = u_view(i, j, k).x() * u_view(i, j, k).x()
                 + u_view(i, j, k).y() * u_view(i, j, k).y()
                 + u_view(i, j, k).z() * u_view(i, j, k).z();


            return weight[vel] * rho_view(i, j, k) * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };

        auto coll = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i )mutable
        {
            if (mesh_view(i.x(), i.y(), i.z()) != 0) {
                for (int vel = 0; vel < Constants->Nvel; vel++) {
                    RealType feq = f_equilibrium(i.x(), i.y(), i.z(), vel);
                    df_post_view(i.x(), i.y(), i.z(), vel) =
                            df_view(i.x(), i.y(), i.z(), vel) - (df_view(i.x(), i.y(), i.z(), vel) - feq) * Constants -> omega;
                }
            }
            else {
                for (int vel = 0; vel < Constants->Nvel; vel++) {
                    df_post_view(i.x(), i.y(), i.z(), vel) = 0;
                }
            }
        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, coll);

    }

    static void streaming(LBMDataPointer &Data, LBMConstantsPointer &Constants) {
        auto mesh_view = Data->meshFluid.getView();
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();



        constexpr int c[27][3] = {
                {0,  0,  0},

                {1,  0,  0},
                {-1, 0,  0},
                {0,  0,  -1},
                {0,  0,  1},
                {0,  -1, 0},
                {0,  1,  0},

                {1,  0,  -1},
                {-1, 0,  1},
                {1,  0,  1},
                {-1, 0,  -1},
                {-1, -1, 0},
                {1,  1,  0},
                {0,  1,  -1},
                {0,  -1, 1},
                {-1, 1,  0},
                {1,  -1, 0},
                {0,  1,  1},
                {0,  -1, -1},

                {-1, 1,  -1},
                {1,  -1, 1},
                {-1, -1, 1},
                {1,  1,  -1},
                {1,  -1, -1},
                {-1, 1,  1},
                {-1, -1, -1},
                {1,  1,  1}};

        auto stream = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i  ) mutable
        {
            if (mesh_view(i.x(), i.y(), i.z()) != 0) { //TODO

                for (int vel = 0; vel < Constants->Nvel; vel++) {
                    int id;
                    int jd;
                    int kd;
                    id = i.x() + c[vel][0];
                    jd = i.y() + c[vel][1];
                    kd = i.z() + c[vel][2];


                    if (id >= 0 && jd >= 0 && kd >= 0 && id < Constants->dimX_int && jd < Constants->dimY_int &&
                        kd < Constants->dimZ_int) {
                        if (mesh_view(id, jd, kd) != 0) {
                            df_view(id, jd, kd, vel) = df_post_view(i.x(), i.y(), i.z(), vel);
                        }
                    }
                }
            }

        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, stream);

    };


    static void bounceBack(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        auto wall_view = Data->meshBoundaryWall;
        auto inlet_view = Data->meshBoundaryInlet;
        auto outlet_view = Data->meshBoundaryOutlet;
        auto mesh_view = Data->meshFluid.getView();
        auto rho_view = Data->rho.getView();
        auto u_view = Data->u.getView();
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();

        constexpr int c[27][3] = {
                {0,  0,  0},

                {1,  0,  0},
                {-1, 0,  0},
                {0,  0,  -1},
                {0,  0,  1},
                {0,  -1, 0},
                {0,  1,  0},

                {1,  0,  -1},
                {-1, 0,  1},
                {1,  0,  1},
                {-1, 0,  -1},
                {-1, -1, 0},
                {1,  1,  0},
                {0,  1,  -1},
                {0,  -1, 1},
                {-1, 1,  0},
                {1,  -1, 0},
                {0,  1,  1},
                {0,  -1, -1},

                {-1, 1,  -1},
                {1,  -1, 1},
                {-1, -1, 1},
                {1,  1,  -1},
                {1,  -1, -1},
                {-1, 1,  1},
                {-1, -1, -1},
                {1,  1,  1}};

        constexpr static RealType weight[27] = {8.f / 27, // 0
                                                2.f / 27, 2.f / 27, 2.f / 27, 2.f / 27, 2.f / 27, 2.f / 27, // 1-6
                                                1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54,
                                                1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, // 7- 18
                                                1.f / 216, 1.f / 216, 1.f / 216, 1.f / 216, 1.f / 216,
                                                1.f / 216, 1.f / 216, 1.f / 216}; // 19-26

        constexpr static int c_rev[27] = {0, 2, 1, 4, 3, 6,
                                          5, 8, 7, 10, 9, 12, 11,
                                          14, 13, 16, 15, 18, 17,
                                          20, 21, 22, 21, 24, 23,
                                          26, 25}; // (reverse velocities are i+1 for i odd and i-1 for i even)

        auto bb_wall = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {


            Vertex vert = wall_view[i.x()].vertex;



            for (int vel = 0; vel < Constants->Nvel; vel++) {

                int dx, dy, dz;

                dx = vert.x - c[vel][0];
                dy = vert.y - c[vel][1];
                dz = vert.z - c[vel][2];

                if (mesh_view(dx, dy, dz) == 0) {
                    df_view(vert.x, vert.y, vert.z, vel) = df_post_view(vert.x, vert.y, vert.z, c_rev[vel]);
                }

            }

        };


        parallelFor<DeviceType>(0, Constants -> wall_num, bb_wall);

        auto bb_inlet = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {
            Vertex vert = inlet_view[i.x()].vertex;
            Vector norm = inlet_view[i.x()].normal;
            Vector velc = inlet_view[i.x()].velocity;

            for (int vel = 0; vel < Constants->Nvel; vel++) {

                if (norm.x() * c[vel][0] + norm.y() * c[vel][1] + norm.z() * c[vel][2] > 0) {
                    df_view(vert.x, vert.y, vert.z, c_rev[vel]) = df_post_view(vert.x, vert.y, vert.z, vel)
                                                                  -
                                                                  6 * weight[vel] * rho_view(vert.x, vert.y, vert.z) * (
                                                                          c[vel][0] * velc.x() + c[vel][1] * velc.y() +
                                                                          c[vel][2] * velc.z()
                                                                  );

                }

            }

        };


        parallelFor<DeviceType>(0, Constants->inlet_num, bb_inlet);

        auto f_equilibrium = [=]
        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k,
        const int &vel ) mutable
        {
            RealType uc, u2;

            uc = c[vel][0] * u_view(i, j, k).x()
                 + c[vel][1] * u_view(i, j, k).y()
                 + c[vel][2] * u_view(i, j, k).z();

            u2 = u_view(i, j, k).x() * u_view(i, j, k).x()
                 + u_view(i, j, k).y() * u_view(i, j, k).y()
                 + u_view(i, j, k).z() * u_view(i, j, k).z();


            return weight[vel] * rho_view(i, j, k) * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };

        auto bb_outlet = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {
            Vertex vert = outlet_view[i.x()].vertex;
            Vector norm = outlet_view[i.x()].normal;
            RealType density = outlet_view[i.x()].density;

            for (int vel = 0; vel < Constants->Nvel; vel++) {

                if (norm.x() * c[vel][0] + norm.y() * c[vel][1] + norm.z() * c[vel][2] > 0) {
                    df_view(vert.x, vert.y, vert.z, c_rev[vel]) = f_equilibrium(vert.x, vert.y, vert.z, vel);

                }

            }

        };


        parallelFor<DeviceType>(0, Constants->outlet_num, bb_outlet);

    }

    static void momentUpdate(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        auto rho_view = Data->rho.getView();
        auto u_view = Data->u.getView();
        auto df_view = Data->df.getView();
        auto mesh_view = Data->meshFluid.getView();


        auto moments = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i )mutable
        {
            if(mesh_view(i.x(), i.y(), i.z()) != 0) {
                rho_view(i.x(), i.y(), i.z()) = df_view(i.x(), i.y(), i.z(), 0) +
                                                df_view(i.x(), i.y(), i.z(), 1) +
                                                df_view(i.x(), i.y(), i.z(), 2) +
                                                df_view(i.x(), i.y(), i.z(), 3) +
                                                df_view(i.x(), i.y(), i.z(), 4) +
                                                df_view(i.x(), i.y(), i.z(), 5) +
                                                df_view(i.x(), i.y(), i.z(), 6) +
                                                df_view(i.x(), i.y(), i.z(), 7) +
                                                df_view(i.x(), i.y(), i.z(), 8) +
                                                df_view(i.x(), i.y(), i.z(), 9) +
                                                df_view(i.x(), i.y(), i.z(), 10) +
                                                df_view(i.x(), i.y(), i.z(), 11) +
                                                df_view(i.x(), i.y(), i.z(), 12) +
                                                df_view(i.x(), i.y(), i.z(), 13) +
                                                df_view(i.x(), i.y(), i.z(), 14) +
                                                df_view(i.x(), i.y(), i.z(), 15) +
                                                df_view(i.x(), i.y(), i.z(), 16) +
                                                df_view(i.x(), i.y(), i.z(), 17) +
                                                df_view(i.x(), i.y(), i.z(), 18) +
                                                df_view(i.x(), i.y(), i.z(), 19) +
                                                df_view(i.x(), i.y(), i.z(), 20) +
                                                df_view(i.x(), i.y(), i.z(), 21) +
                                                df_view(i.x(), i.y(), i.z(), 22) +
                                                df_view(i.x(), i.y(), i.z(), 23) +
                                                df_view(i.x(), i.y(), i.z(), 24) +
                                                df_view(i.x(), i.y(), i.z(), 25) +
                                                df_view(i.x(), i.y(), i.z(), 26);

                u_view(i.x(), i.y(), i.z())(0) = (df_view(i.x(), i.y(), i.z(), 1)
                                                  - df_view(i.x(), i.y(), i.z(), 2)
                                                  + df_view(i.x(), i.y(), i.z(), 7)
                                                  - df_view(i.x(), i.y(), i.z(), 8)
                                                  + df_view(i.x(), i.y(), i.z(), 9)
                                                  - df_view(i.x(), i.y(), i.z(), 10)
                                                  - df_view(i.x(), i.y(), i.z(), 11)
                                                  + df_view(i.x(), i.y(), i.z(), 12)
                                                  - df_view(i.x(), i.y(), i.z(), 15)
                                                  + df_view(i.x(), i.y(), i.z(), 16)
                                                  - df_view(i.x(), i.y(), i.z(), 19)
                                                  + df_view(i.x(), i.y(), i.z(), 20)
                                                  - df_view(i.x(), i.y(), i.z(), 21)
                                                  + df_view(i.x(), i.y(), i.z(), 22)
                                                  + df_view(i.x(), i.y(), i.z(), 23)
                                                  - df_view(i.x(), i.y(), i.z(), 24)
                                                  - df_view(i.x(), i.y(), i.z(), 25)
                                                  + df_view(i.x(), i.y(), i.z(), 26)) / rho_view(i.x(), i.y(), i.z());


                u_view(i.x(), i.y(), i.z())(1) = (-df_view(i.x(), i.y(), i.z(), 5)
                                                  + df_view(i.x(), i.y(), i.z(), 6)
                                                  - df_view(i.x(), i.y(), i.z(), 11)
                                                  + df_view(i.x(), i.y(), i.z(), 12)
                                                  + df_view(i.x(), i.y(), i.z(), 13)
                                                  - df_view(i.x(), i.y(), i.z(), 14)
                                                  + df_view(i.x(), i.y(), i.z(), 15)
                                                  - df_view(i.x(), i.y(), i.z(), 16)
                                                  + df_view(i.x(), i.y(), i.z(), 17)
                                                  - df_view(i.x(), i.y(), i.z(), 18)
                                                  + df_view(i.x(), i.y(), i.z(), 19)
                                                  - df_view(i.x(), i.y(), i.z(), 20)
                                                  - df_view(i.x(), i.y(), i.z(), 21)
                                                  + df_view(i.x(), i.y(), i.z(), 22)
                                                  - df_view(i.x(), i.y(), i.z(), 23)
                                                  + df_view(i.x(), i.y(), i.z(), 24)
                                                  - df_view(i.x(), i.y(), i.z(), 25)
                                                  + df_view(i.x(), i.y(), i.z(), 26)) / rho_view(i.x(), i.y(), i.z());

                u_view(i.x(), i.y(), i.z())(2) = (-df_view(i.x(), i.y(), i.z(), 3)
                                                  + df_view(i.x(), i.y(), i.z(), 4)
                                                  - df_view(i.x(), i.y(), i.z(), 7)
                                                  + df_view(i.x(), i.y(), i.z(), 8)
                                                  + df_view(i.x(), i.y(), i.z(), 9)
                                                  - df_view(i.x(), i.y(), i.z(), 10)
                                                  - df_view(i.x(), i.y(), i.z(), 13)
                                                  + df_view(i.x(), i.y(), i.z(), 14)
                                                  + df_view(i.x(), i.y(), i.z(), 17)
                                                  - df_view(i.x(), i.y(), i.z(), 18)
                                                  - df_view(i.x(), i.y(), i.z(), 19)
                                                  + df_view(i.x(), i.y(), i.z(), 20)
                                                  + df_view(i.x(), i.y(), i.z(), 21)
                                                  - df_view(i.x(), i.y(), i.z(), 22)
                                                  - df_view(i.x(), i.y(), i.z(), 23)
                                                  + df_view(i.x(), i.y(), i.z(), 24)
                                                  - df_view(i.x(), i.y(), i.z(), 25)
                                                  + df_view(i.x(), i.y(), i.z(), 26)) / rho_view(i.x(), i.y(), i.z());
            }
        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, moments);


    }

    static void errorEvaluation(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        auto uStatArr_view = Data->u.getStorageArray().getConstView();
        auto uErrorStatArr_view = Data->u_error.getStorageArray().getConstView();
        auto mesh_view = Data->meshFluid.getStorageArray().getConstView();

        auto u_view = Data->u.getView();
        auto uError_view = Data->u_error.getView();
        

        
        RealType err1, err2;
        err1 = 0.0f;
        err2 = 0.0f;

        auto fetch1 = [=]
        __cuda_callable__(int
        i ) ->RealType
        {
            //dui^2
            if (mesh_view[i] == 1 ){
                const RealType &add1 = (uStatArr_view[i].x() - uErrorStatArr_view[i].x()) * (uStatArr_view[i].x() - uErrorStatArr_view[i].x()) +
                                       (uStatArr_view[i].y() - uErrorStatArr_view[i].y()) * (uStatArr_view[i].y() - uErrorStatArr_view[i].y()) +
                                       (uStatArr_view[i].y() - uErrorStatArr_view[i].z()) * (uStatArr_view[i].y() - uErrorStatArr_view[i].z());

                return sqrt(add1);
            }
            return 0;
        };
        
        auto fetch2 = [=]
        __cuda_callable__(int
        i ) ->RealType
        {   
            //ui^2
            if (mesh_view[i] == 1) {
                const RealType &add2 =
                        uStatArr_view[i].x() * uStatArr_view[i].x() +
                        uStatArr_view[i].y() * uStatArr_view[i].y() +
                        uStatArr_view[i].z() * uStatArr_view[i].z();

                return sqrt(add2);
            }
            return 0;
        };
        
        err1 = sqrt(reduce<DeviceType, int>(0,
                                            uStatArr_view.getSize(),
                                            fetch1,
                                            TNL::Plus(),
                                            0.0));

        err2 = sqrt(reduce<DeviceType, int>(0,
                                            uStatArr_view.getSize(),
                                            fetch2,
                                            TNL::Plus(),
                                            0.0));
        Constants -> err =err1/err2;

        auto copy = [=] __cuda_callable__ ( const TNL::Containers::StaticArray< 3, int >& i  ) mutable
        {
            uError_view(i.x(), i.y(), i.z()) = u_view(i.x(), i.y(), i.z());

        };

        TNL::Containers::StaticArray< 3, int > begin{ 0, 0, 0};
        TNL::Containers::StaticArray< 3, int > end{ Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor< DeviceType >( begin, end, copy );
    }

    auto velocityInletAverage(LBMDataPointer &Data, LBMConstantsPointer &Constants) {
/*todo
        auto inlet_view = Data -> meshBoundaryInlet.getView();

        auto fetch_u = [=]
        __cuda_callable__(
        const int &i)
        {
            return inlet_view(i).velocity;

        };

        /*auto reduce_u = [=](
                const VectorType& a, VectorType& b){

            return a + b;
        };

        VectorType u_avrg = reduce<DeviceType, int>(0,
                                                    Constants -> inlet_num,
                                                    fetch_u,
                                                    TNL::Plus,
                                                    0.f );
        return u_avrg;
*/
    };

    auto densityOutletAverage() {
        /* todo
            /*auto fetch_rho = [=]
           __cuda_callable__(
           const int &i)
           {
               return rho_view(i).density ;
           };

           auto reduce_rho = [=](
                   const RealType &a, const RealType &b){


               return a+b;
           };*/

    };

    const static int numberOfDiscreteVelocities = 27;
    enum {
        //velocities are denoted by 4 cardinal directions:
        // "n" - north, "s" - south, "e" - east and "w" - west
        // in addition with "f" - forward, "b" - backward direction,
        // "o" - non-part of vector.
        //
        //           +N2
        //            ↑   -B3
        //            | ⋰
        //  -W1 <-----O-----> +E1
        //          ⋰ |
        //      +F3   ↓
        //           -S2
        //
        //            +y
        //            ↑    -z
        //            | ⋰
        //   -x <-----O-----> +x
        //          ⋰ |
        //      +z    ↓
        //            -y
        //
        // index precedence is:
        // 1: {w, o, e}
        // 2: {s, o, n}
        // 3: {b, o, f}

        // weight is 8/27 for 0
        ooo = 0,

        // weight is 2/27 for (1-6)
        eoo = 1,
        woo = 2,
        oob = 3,
        oof = 4,
        oso = 5,
        ono = 6,

        // weight is 1/54 for (7-18)
        eob = 7,
        wof = 8,
        eof = 9,
        wob = 10,
        wso = 11,
        eno = 12,
        onb = 13,
        osf = 14,
        wno = 15,
        eso = 16,
        onf = 17,
        osb = 18,

        // weight is 1/126 for (19-26)
        wnb = 19,
        esf = 20,
        wsf = 21,
        enb = 22,
        esb = 23,
        wnf = 24,
        wsb = 25,
        enf = 26

    };

//const static int Nvel = 27; //number of discrete velocities
    constexpr static int c[27][3] = {
            {0,  0,  0},

            {1,  0,  0},
            {-1, 0,  0},
            {0,  0,  -1},
            {0,  0,  1},
            {0,  -1, 0},
            {0,  1,  0},

            {1,  0,  -1},
            {-1, 0,  1},
            {1,  0,  1},
            {-1, 0,  -1},
            {-1, -1, 0},
            {1,  1,  0},
            {0,  1,  -1},
            {0,  -1, 1},
            {-1, 1,  0},
            {1,  -1, 0},
            {0,  1,  1},
            {0,  -1, -1},

            {-1, 1,  -1},
            {1,  -1, 1},
            {-1, -1, 1},
            {1,  1,  -1},
            {1,  -1, -1},
            {-1, 1,  1},
            {-1, -1, -1},
            {1,  1,  1}};

    constexpr static RealType weight[27] = {8.f / 27, // 0
                                            2.f / 27, 2.f / 27, 2.f / 27, 2.f / 27, 2.f / 27, 2.f / 27, // 1-6
                                            1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54,
                                            1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, 1.f / 54, // 7- 18
                                            1.f / 216, 1.f / 216, 1.f / 216, 1.f / 216, 1.f / 216,
                                            1.f / 216, 1.f / 216, 1.f / 216}; // 19-26

    constexpr static int c_rev[27] = {0, 2, 1, 4, 3, 6,
                                      5, 8, 7, 10, 9, 12, 11,
                                      14, 13, 16, 15, 18, 17,
                                      20, 21, 22, 21, 24, 23,
                                      26, 25}; // (reverse velocities are i+1 for i odd and i-1 for i even)


};

#endif //D3Q27_H
