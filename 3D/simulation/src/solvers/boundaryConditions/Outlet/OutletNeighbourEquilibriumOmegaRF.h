//
// Created by stloufra on 10/30/23.
//

#ifndef OUTLETNEIGHBOUREQUILIBRIUMOMEGARF_H
#define OUTLETNEIGHBOUREQUILIBRIUMOMEGARF_H

#include "../../../traits/LBMTraits.h"

template<typename MODELDATA>
struct OutletNeighbourEquilibriumOmegaRF{
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;



    static void outletOmega(LBMDataPointer &Data, LBMConstantsPointer &Constants) {
        // Dumping by omega 10% of the length in normal direction in front of the outlet
        // by addition of linear function from omegaDumpingLow to omegaDumpingHigh
        // (stored in Constants and needs to be intialized in main).
        // So far implemented only for outlets at ends of the domain in
        // any direction tangential to one of the axis.

        auto outletHost_view = Data->meshBoundaryOutletHost.getView();

        Vector norm = outletHost_view[0].normal;

        Vector normPlane;

        normPlane.x() = abs(abs(norm.x()) - 1);
        normPlane.y() = abs(abs(norm.y()) - 1);
        normPlane.z() = abs(abs(norm.z()) - 1); // normal plane perpendicular to the normal vector


        int dimX_01_int = std::ceil(Constants->dimX_int * 0.1f);
        int dimY_01_int = std::ceil(Constants->dimY_int * 0.1f);
        int dimZ_01_int = std::ceil(Constants->dimZ_int * 0.1f); // 10% of lenghts

        int dimX = dimX_01_int * abs(norm.x()) + Constants->dimX_int * abs(normPlane.x());
        int dimY = dimY_01_int * abs(norm.y()) + Constants->dimY_int * abs(normPlane.y());
        int dimZ = dimZ_01_int * abs(norm.z()) + Constants->dimZ_int * abs(normPlane.z());
        // dimensions of domain to be dumped according to normal axis

        int beginX = Constants->dimX_int - dimX;
        int beginY = Constants->dimY_int - dimY;
        int beginZ = Constants->dimZ_int - dimZ; // possible beginnings

        int intervalSize = abs(dimX_01_int * norm.x() + dimY_01_int * norm.y() + dimZ_01_int * norm.z());

        RealType omegaDumpingLow =  Constants->omegaDumpingLow;
        RealType omegaDumpingHigh =  Constants->omegaDumpingHigh;

        RealType derivation = (omegaDumpingHigh - omegaDumpingLow) / intervalSize;


        ///////////////////////////////////////////////////////////////////////////////////////////

        auto omega_view = Data->omega.getView();


        auto omegaFunc = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i ) mutable
        {
            int possition =     (i.x() - beginX) * norm.x()* norm.x() +
                                (i.y() - beginY) * norm.y()* norm.y() +
                                (i.z() - beginZ) * norm.z()* norm.z();

            RealType omegaDumping = derivation * possition;

            omega_view(i.x(), i.y(), i.z()) = omega_view(i.x(), i.y(), i.z()) - omegaDumping;

        };
        TNL::Containers::StaticArray<3, int> begin1{beginX, beginY, beginZ};
        TNL::Containers::StaticArray<3, int> end1{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin1, end1, omegaFunc);

    }

    static void outlet(LBMDataPointer &Data, LBMConstantsPointer &Constants) {

        // f = feq + fneqNeigh
        // when u_x < 0 then p = p_s; u_x > 0 then p = p_s - 0.5*rho*u^2


        auto outlet_view = Data->meshBoundaryOutlet.getView();
        auto u_view = Data->u.getView();

        auto rho_view = Data->rho.getView();
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();
        auto mesh_view = Data->meshFluid.getView();

        const auto Nvel = Constants->Nvel;

        MODELDATA MD;

        auto f_equilibrium_outlet = [=]
        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k,
        const int &vel,
        const RealType &density) mutable
        {
            RealType uc, u2;

            uc = MD.c[vel][0] * u_view(i, j, k).x()
                 + MD.c[vel][1] * u_view(i, j, k).y()
                 + MD.c[vel][2] * u_view(i, j, k).z();

            u2 = u_view(i, j, k).x() * u_view(i, j, k).x()
                 + u_view(i, j, k).y() * u_view(i, j, k).y()
                 + u_view(i, j, k).z() * u_view(i, j, k).z();


            return MD.weight[vel] * density * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };

        auto f_equilibrium = [=]
        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k,
        const int &vel ) mutable
        {
            RealType uc, u2;

            uc = MD.c[vel][0] * u_view(i, j, k).x()
                 + MD.c[vel][1] * u_view(i, j, k).y()
                 + MD.c[vel][2] * u_view(i, j, k).z();

            u2 = u_view(i, j, k).x() * u_view(i, j, k).x()
                 + u_view(i, j, k).y() * u_view(i, j, k).y()
                 + u_view(i, j, k).z() * u_view(i, j, k).z();


            return MD.weight[vel] * rho_view(i, j, k) * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };

        auto bb_outlet = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &i  ) mutable
        {
            Vertex vert = outlet_view[i.x()].vertex;
            Vector norm = outlet_view[i.x()].normal;
            RealType density = outlet_view[i.x()].density;

            if(u_view(vert.x, vert.y, vert.z).x() < 0 )
            {
                density = density - 0.5*u_view(vert.x, vert.y, vert.z).x()*u_view(vert.x, vert.y, vert.z).x()*density;
            }


            bool reg = outlet_view[i.x()].regular;

            for (int vel = 0; vel < Nvel; vel++) {


                if (norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0) {

                    if (reg) {


                        int xNeigh = vert.x + MD.c[vel][0];
                        int yNeigh = vert.y + MD.c[vel][1];
                        int zNeigh = vert.z + MD.c[vel][2];

                        RealType fEqNeigh = f_equilibrium(xNeigh, yNeigh, zNeigh, vel);
                        RealType fNeqNeigh = df_view(xNeigh, yNeigh, zNeigh, vel);
                        RealType fEqBC = f_equilibrium_outlet(vert.x, vert.y, vert.z, vel, density);

                        df_view(vert.x, vert.y, vert.z, vel) = fEqBC + (fEqNeigh - fNeqNeigh);
                    } else {
                        int xNeigh = vert.x + MD.c[vel][0];
                        int yNeigh = vert.y + MD.c[vel][1];
                        int zNeigh = vert.z + MD.c[vel][2];

                        if (mesh_view(xNeigh, yNeigh, zNeigh) > 0) {

                            RealType fEqNeigh = f_equilibrium(xNeigh, yNeigh, zNeigh, vel);
                            RealType fNeqNeigh = df_view(xNeigh, yNeigh, zNeigh, vel);
                            RealType fEqBC = f_equilibrium_outlet(vert.x, vert.y, vert.z, vel, density);

                            df_view(vert.x, vert.y, vert.z, vel) = fEqBC + (fEqNeigh - fNeqNeigh);
                        }

                    }
                }
            }
        };


        parallelFor<DeviceType>(0, Constants->outlet_num, bb_outlet);

    }

};

#endif
