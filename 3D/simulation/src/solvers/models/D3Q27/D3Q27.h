#ifndef D3Q27_H
#define  D3Q27_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>
#include <omp.h>
#include "../../../geometry/geometryMesherBoundary.h"
#include "../../../traits/LBMTraits.h"
#include "D3Q27Data.h"
#include "D3Q27CollisionRolled.h"
#include "D3Q27CollisionUnrolled.h"
#include "D3Q27InitVectorEq.h"
#include "D3Q27InitVariablesEq.h"
#include <TNL/Containers/StaticArray.h>
#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Algorithms/reduce.h>

using namespace TNL;
using namespace TNL::Algorithms;

template<typename MODELDATATYPE, typename COLLISIONTYPE, typename INITIALIZATIONTYPE> //TODO implement typename COL
class D3Q27 {
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;


public:

    static void initialization(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        INITIALIZATIONTYPE::initialization(Data, Constants);
    };

    static void collision(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        COLLISIONTYPE::collision(Data, Constants);

    }


    static void streaming(LBMDataPointer &Data, LBMConstantsPointer &Constants) {
        auto mesh_view = Data->meshFluid.getView();
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();

        const auto Nvel = MODELDATATYPE::numberOfDiscreteVelocities;
        auto dimX_int = Constants->dimX_int;
        auto dimY_int = Constants->dimY_int;
        auto dimZ_int = Constants->dimZ_int;

        D3Q27Data MDATA;

        auto stream = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i  ) mutable
        {
            if (mesh_view(i.x(), i.y(), i.z()) != 0) { //TODO

                for (int vel = 0; vel < Nvel; vel++) {
                    int id;
                    int jd;
                    int kd;
                    id = i.x() + MDATA.c[vel][0];
                    jd = i.y() + MDATA.c[vel][1];
                    kd = i.z() + MDATA.c[vel][2];


                    if (id >= 0 && jd >= 0 && kd >= 0 && id < dimX_int && jd < dimY_int &&
                        kd < dimZ_int) {
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

        const auto Nvel = MODELDATATYPE::numberOfDiscreteVelocities;

        D3Q27Data MDATA;


        auto bb_wall = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {


            Vertex vert = wall_view[i.x()].vertex;


            for (int vel = 0; vel < Nvel; vel++) {

                int dx, dy, dz;

                dx = vert.x - MDATA.c[vel][0];
                dy = vert.y - MDATA.c[vel][1];
                dz = vert.z - MDATA.c[vel][2];

                if (mesh_view(dx, dy, dz) == 0) {
                    df_view(vert.x, vert.y, vert.z, vel) = df_post_view(vert.x, vert.y, vert.z, MDATA.c_rev[vel]);
                }

            }

        };


        parallelFor<DeviceType>(0, Constants->wall_num, bb_wall);

        auto bb_inlet = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {
            Vertex vert = inlet_view[i.x()].vertex;
            Vector norm = inlet_view[i.x()].normal;
            Vector velc = inlet_view[i.x()].velocity;


            for (int vel = 0; vel < Nvel; vel++) {

                if (norm.x() * MDATA.c[vel][0] + norm.y() * MDATA.c[vel][1] + norm.z() * MDATA.c[vel][2] > 0) {
                    df_view(vert.x, vert.y, vert.z, MDATA.c_rev[vel]) = df_post_view(vert.x, vert.y, vert.z, vel)
                                                                        - 6 * MDATA.weight[vel] *
                                                                          rho_view(vert.x, vert.y, vert.z) * (
                                                                                  MDATA.c[vel][0] * velc.x() +
                                                                                  MDATA.c[vel][1] * velc.y() +
                                                                                  MDATA.c[vel][2] * velc.z()
                                                                          );

                }

            }

        };


        parallelFor<DeviceType>(0, Constants->inlet_num, bb_inlet);

        auto f_equilibrium_outlet = [=]
        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k,
        const int &vel,
        const RealType &density) mutable
        {
            RealType uc, u2;

            uc = MDATA.c[vel][0] * u_view(i, j, k).x()
                 + MDATA.c[vel][1] * u_view(i, j, k).y()
                 + MDATA.c[vel][2] * u_view(i, j, k).z();

            u2 = u_view(i, j, k).x() * u_view(i, j, k).x()
                 + u_view(i, j, k).y() * u_view(i, j, k).y()
                 + u_view(i, j, k).z() * u_view(i, j, k).z();


            return MDATA.weight[vel] * density * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };

        auto bb_outlet = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {
            Vertex vert = outlet_view[i.x()].vertex;
            Vector norm = outlet_view[i.x()].normal;
            RealType density = outlet_view[i.x()].density;

            for (int vel = 0; vel < Nvel; vel++) {

                if (norm.x() * MDATA.c[vel][0] + norm.y() * MDATA.c[vel][1] + norm.z() * MDATA.c[vel][2] > 0) {
                    df_view(vert.x, vert.y, vert.z, MDATA.c_rev[vel]) = f_equilibrium_outlet(vert.x, vert.y, vert.z,vel, density );

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
            if (mesh_view(i.x(), i.y(), i.z()) != 0) {
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


    static void pressureUpdate(LBMDataPointer &Data, LBMConstantsPointer &Constants) {


        auto p_view = Data->p.getView();
        auto rho_view = Data->rho.getView();
        auto mesh_view = Data->meshFluid.getView();
        auto cs2 = Constants-> cs2;


        auto pressure = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i)
        mutable
        {

            p_view(i.x(), i.y(), i.z()) = rho_view(i.x(), i.y(), i.z())*cs2;

        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, pressure);
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
            if (mesh_view[i] == 1) {
                const RealType &add1 = (uStatArr_view[i].x() - uErrorStatArr_view[i].x()) *
                                       (uStatArr_view[i].x() - uErrorStatArr_view[i].x()) +
                                       (uStatArr_view[i].y() - uErrorStatArr_view[i].y()) *
                                       (uStatArr_view[i].y() - uErrorStatArr_view[i].y()) +
                                       (uStatArr_view[i].y() - uErrorStatArr_view[i].z()) *
                                       (uStatArr_view[i].y() - uErrorStatArr_view[i].z());

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
        Constants->err = err1 / err2;

        auto copy = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i  ) mutable
        {
            uError_view(i.x(), i.y(), i.z()) = u_view(i.x(), i.y(), i.z());

        };

        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, copy);
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
           };
        */
    };

};

#endif //D3Q27_H
