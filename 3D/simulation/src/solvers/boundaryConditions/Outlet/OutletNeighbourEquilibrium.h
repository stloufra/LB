//
// Created by stloufra on 10/30/23.
//

#ifndef OUTLETNEIGHBOUREQUILIBRIUM_H
#define OUTLETNEIGHBOUREQUILIBRIUM_H

#include "../../../traits/LBMTraits.h"

template<typename MODELDATA>
struct OutletNeighbourEquilibrium {
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void outletOmega(LBMDataPointer &Data, LBMConstantsPointer &Constants) {}

    static void outlet(LBMDataPointer &Data, LBMConstantsPointer &Constants) {


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
            bool reg = outlet_view[i.x()].regular;

            for (int vel = 0; vel < Nvel; vel++) {


                if(norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0) {

                    if(reg) {


                        int xNeigh = vert.x + MD.c[vel][0];
                        int yNeigh = vert.y + MD.c[vel][1];
                        int zNeigh = vert.z + MD.c[vel][2];

                        RealType fEqNeigh = f_equilibrium(xNeigh, yNeigh, zNeigh, vel);
                        RealType fNeqNeigh = df_view(xNeigh, yNeigh, zNeigh, vel);
                        RealType fEqBC = f_equilibrium_outlet(vert.x, vert.y, vert.z, vel, density);

                        df_view(vert.x, vert.y, vert.z, vel) = fEqBC + (fEqNeigh - fNeqNeigh);
                    }
                    else{
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
