//
// Created by stloufra on 10/30/23.
//

#ifndef OUTLETDENSITYINTERPOLATED_H
#define OUTLETDENSITYINTERPOLATED_H

#include "../../../traits/LBMTraits.h"

template<typename MODELDATA>
struct OutletDensityInterpolated {
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void outletOmega(LBMDataPointer &Data, LBMConstantsPointer &Constants) {}

    static void outlet(LBMDataPointer &Data, LBMConstantsPointer &Constants) {


        auto outlet_view = Data->meshBoundaryOutlet.getView();
        auto mesh_view = Data->meshFluid.getView();
        auto u_view = Data->u.getView();

        auto rho_view = Data->rho.getView();
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();

        const auto Nvel = Constants->Nvel;
        const auto cs = Constants->cs;


        MODELDATA MD;

        auto f_equilibrium_defined = [=]
        __cuda_callable__(
        const RealType &ux,
        const RealType &uy,
        const RealType &uz,
        const RealType &density,
        const int &vel) mutable
        {
            RealType uc, u2;

            uc = MD.c[vel][0] *ux + MD.c[vel][1] *uy + MD.c[vel][2] *uz;

            u2 =ux *ux +uy *uy +uz *uz;

            return MD.weight[vel] * density * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };
        
        auto computeRho = [=]
        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k) mutable
        {
          return     df_view(i, j, k, 0) +
                    df_view(i, j, k, 1) +
                    df_view(i, j, k, 2) +
                    df_view(i, j, k, 3) +
                    df_view(i, j, k, 4) +
                    df_view(i, j, k, 5) +
                    df_view(i, j, k, 6) +
                    df_view(i, j, k, 7) +
                    df_view(i, j, k, 8) +
                    df_view(i, j, k, 9) +
                    df_view(i, j, k, 10) +
                    df_view(i, j, k, 11) +
                    df_view(i, j, k, 12) +
                    df_view(i, j, k, 13) +
                    df_view(i, j, k, 14) +
                    df_view(i, j, k, 15) +
                    df_view(i, j, k, 16) +
                    df_view(i, j, k, 17) +
                    df_view(i, j, k, 18) +
                    df_view(i, j, k, 19) +
                    df_view(i, j, k, 20) +
                    df_view(i, j, k, 21) +
                    df_view(i, j, k, 22) +
                    df_view(i, j, k, 23) +
                    df_view(i, j, k, 24) +
                    df_view(i, j, k, 25) +
                    df_view(i, j, k, 26);
          };

        auto computeUx = [=]
        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k,
        const RealType &density) mutable
        {
            return     (df_view(i, j, k, 1)
                      - df_view(i, j, k, 2)
                      + df_view(i, j, k, 7)
                      - df_view(i, j, k, 8)
                      + df_view(i, j, k, 9)
                      - df_view(i, j, k, 10)
                      - df_view(i, j, k, 11)
                      + df_view(i, j, k, 12)
                      - df_view(i, j, k, 15)
                      + df_view(i, j, k, 16)
                      - df_view(i, j, k, 19)
                      + df_view(i, j, k, 20)
                      - df_view(i, j, k, 21)
                      + df_view(i, j, k, 22)
                      + df_view(i, j, k, 23)
                      - df_view(i, j, k, 24)
                      - df_view(i, j, k, 25)
                      + df_view(i, j, k, 26)) / density;
        };

        auto computeUy = [=]
        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k,
        const RealType &density) mutable
        {
            return    (-df_view(i, j, k, 5)
                      + df_view(i, j, k, 6)
                      - df_view(i, j, k, 11)
                      + df_view(i, j, k, 12)
                      + df_view(i, j, k, 13)
                      - df_view(i, j, k, 14)
                      + df_view(i, j, k, 15)
                      - df_view(i, j, k, 16)
                      + df_view(i, j, k, 17)
                      - df_view(i, j, k, 18)
                      + df_view(i, j, k, 19)
                      - df_view(i, j, k, 20)
                      - df_view(i, j, k, 21)
                      + df_view(i, j, k, 22)
                      - df_view(i, j, k, 23)
                      + df_view(i, j, k, 24)
                      - df_view(i, j, k, 25)
                      + df_view(i, j, k, 26)) / density;
        };

        auto computeUz = [=]
        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k,
        const RealType &density) mutable
        {
            return   (-df_view(i, j, k, 3)
                      + df_view(i, j, k, 4)
                      - df_view(i, j, k, 7)
                      + df_view(i, j, k, 8)
                      + df_view(i, j, k, 9)
                      - df_view(i, j, k, 10)
                      - df_view(i, j, k, 13)
                      + df_view(i, j, k, 14)
                      + df_view(i, j, k, 17)
                      - df_view(i, j, k, 18)
                      - df_view(i, j, k, 19)
                      + df_view(i, j, k, 20)
                      + df_view(i, j, k, 21)
                      - df_view(i, j, k, 22)
                      - df_view(i, j, k, 23)
                      + df_view(i, j, k, 24)
                      - df_view(i, j, k, 25)
                      + df_view(i, j, k, 26)) / density;
        };

        auto bb_outlet = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &nod  ) mutable
        {
            Vertex vert = outlet_view[nod.x()].vertex;
            Vector norm = outlet_view[nod.x()].normal;
            RealType density = outlet_view[nod.x()].density;
            bool reg = outlet_view[nod.x()].regular;

            int i = vert.x;
            int j = vert.y;
            int k = vert.z;

            Vector xp(1.f, 0.f, 0.f), xm(-1.f, 0.f, 0.f);

            if ( norm == xp){
                if(!reg)
                {
                    for (int vel = 0; vel < Nvel; vel++){

                        int dy = vert.y - MD.c[vel][1];
                        int dz = vert.z - MD.c[vel][2];

                        if(mesh_view(i, dy, dz) == 0)
                        {
                            df_view(i, j, k, vel) = df_view(i, j, k, MD.c_rev[vel]);
                        }
                    }
                }

                for (int vel = 0; vel < Nvel; vel++) {
                    if (norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0) { //pokud miri dovnitr

                        int dx;

                        dx = i - 1;

                        df_view(i,j,k,vel)=cs*df_post_view(dx,j,k,vel) + (1-cs)*df_post_view(i,j,k, vel);

                    }
                }

                RealType rho = computeRho(i,j,k);
                RealType ux = computeUx(i,j,k,density);
                RealType uy = computeUy(i,j,k,density);
                RealType uz = computeUz(i,j,k,density);

                for (int vel = 0; vel < Nvel; vel++) {

                    df_view(i,j,k,vel) = df_view(i,j,k,vel) - f_equilibrium_defined(ux, uy, uz, rho, vel) + f_equilibrium_defined(ux, uy, uz, 1.f, vel);

                }
            }
            else if ( norm == xm){
                for (int vel = 0; vel < Nvel; vel++) {
                    if (norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0) { //pokud miri dovnitr

                        int dx;

                        dx = i + 1;

                        df_view(i,j,k,vel)=cs*df_post_view(dx,j,k,vel) + (1-cs)*df_post_view(i,j,k, vel);

                    }
                }

                RealType rho = computeRho(i,j,k);
                RealType ux = computeUx(i,j,k,density);
                RealType uy = computeUy(i,j,k,density);
                RealType uz = computeUz(i,j,k,density);

                for (int vel = 0; vel < Nvel; vel++) {

                    df_view(i,j,k,vel) = df_view(i,j,k,vel) - f_equilibrium_defined(ux, uy, uz, rho, vel) + f_equilibrium_defined(ux, uy, uz, 1.f, vel);

                }

            }
            else{
                printf("Not yet supported outlet.\n");
            }

        };


        parallelFor<DeviceType>(0, Constants->outlet_num, bb_outlet);

    }

};

#endif
