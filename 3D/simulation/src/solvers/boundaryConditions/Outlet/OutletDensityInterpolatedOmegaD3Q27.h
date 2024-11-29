//
// Created by stloufra on 10/30/23.
//

#ifndef OUTLETDENSITYINTERPOLATEDOMEGAD3Q27_H
#define OUTLETDENSITYINTERPOLATEDOMEGAD3Q27_H

#include "../../../traits/LBMTraits.h"

template <typename MODELDATA>
struct OutletDensityInterpolatedOmegaD3Q27
{
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


        int dimX_01_int = std::ceil(Constants->dimX_int * 0.05f);
        int dimY_01_int = std::ceil(Constants->dimY_int * 0.05f);
        int dimZ_01_int = std::ceil(Constants->dimZ_int * 0.05f); // 5% of lenghts

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


    static void outlet(LBMDataPointer& Data, LBMConstantsPointer& Constants)
    {
        auto outlet_view = Data->meshBoundaryOutlet.getView();
        auto mesh_view = Data->meshFluid.getView();
        auto ux_view = Data->ux.getView();

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

        const int& vel
        )
        mutable
        {
            RealType uc, u2;

            uc = MD.c[vel][0] * ux + MD.c[vel][1] * uy + MD.c[vel][2] * uz;

            u2 = ux * ux + uy * uy + uz * uz;

            return MD.weight[vel] * density * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };

        auto computeRho = [=]

        __cuda_callable__(
        const int &i,
        const int &j,
        const int &k)
        mutable
        {
            return df_view(i, j, k, 0) +
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

        const RealType& density
        )
        mutable
        {
            return (df_view(i, j, k, 1)
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

        const RealType& density
        )
        mutable
        {
            return (-df_view(i, j, k, 5)
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

        const RealType& density
        )
        mutable
        {
            return (-df_view(i, j, k, 3)
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

        const TNL::Containers::StaticArray<1, int>& nod
        )
        mutable
        {
            Vertex vert = outlet_view[nod.x()].vertex;
            Vector norm = outlet_view[nod.x()].normal;
            RealType density = outlet_view[nod.x()].density;
            bool reg = outlet_view[nod.x()].regular;

            int i = vert.x;
            int j = vert.y;
            int k = vert.z;

            Vector xp(1.f, 0.f, 0.f), xm(-1.f, 0.f, 0.f);

            if (norm == xp)
            {
                if (!reg)
                {
                    for (int vel = 0; vel < Nvel; vel++) //BBstěna
                    {
                        int dy = vert.y - MD.c[vel][1];
                        int dz = vert.z - MD.c[vel][2];

                        if (mesh_view(i, dy, dz) == 0)
                        {
                            df_view(i, j, k, vel) = df_post_view(i, j, k, MD.c_rev[vel]);
                        }
                    }

                    /*if((mesh_view(i-1, j, k) == -3)) //proti sm2ru normaly je periodicita
                    {
                        for (int vel = 0; vel < Nvel; vel++)
                        {
                            int dx = vert.y - MD.c[vel][0];
                            int dy = vert.y - MD.c[vel][1];
                            int dz = vert.z - MD.c[vel][2];

                            if(mesh_view(dx, dy, dz) == 0) // a rychlost zatim nevim od kud predepsat
                            {
                                if(j=16)
                                {
                                    df_view(i, j, k, vel) = df_post_view(i, 1, k, vel);
                                }else if(j=1)
                                {
                                    df_view(i, j, k, vel) = df_post_view(i, 16, k, vel);
                                }
                            }
                        }
                    }*/

                    for (int vel = 0; vel < Nvel; vel++)
                    {
                        if (j == 16)
                        {
                            if (MD.c[vel][1] < 0) //periodicita
                            {
                                df_view(i, j, k, vel) = df_view(i, 0, k, vel);
                            }

                            if (k == 220 && MD.c[vel][2] < 0) //odraz nahore
                            {
                                if(MD.c[vel][1] <= 0){
                                    df_view(i, j, k, vel) = df_post_view(i, j, k, MD.c_rev[vel]);
                                }else if(MD.c[vel][1] > 0){
                                    df_view(i, j, k, vel) = df_post_view(i, 1, k, MD.c_rev[vel]);
                                }
                            }
                            else if (k == 1 && MD.c[vel][2] > 0) //odraz dole
                            {
                                if(MD.c[vel][1] <= 0){
                                    df_view(i, j, k, vel) = df_post_view(i, j, k, MD.c_rev[vel]);
                                }else if(MD.c[vel][1] > 0){
                                    df_view(i, j, k, vel) = df_post_view(i, 1, k, MD.c_rev[vel]);
                                }
                            }
                        }
                        if (j == 1)
                        {
                            if (MD.c[vel][1] > 0) //periodicita
                            {
                                df_view(i, j, k, vel) = df_view(i, 17, k, vel);
                            }

                            if (k == 220 && MD.c[vel][2] < 0)//odraz nahore
                            {
                                if(MD.c[vel][1] >= 0){
                                    df_view(i, j, k, vel) = df_post_view(i, j, k, MD.c_rev[vel]);
                                }else if(MD.c[vel][1] < 0){
                                    df_view(i, j, k, vel) = df_post_view(i, 16, k, MD.c_rev[vel]);
                                }
                            }
                            else if (k == 1 && MD.c[vel][2] > 0)
                            {
                                if(MD.c[vel][1] >= 0){
                                    df_view(i, j, k, vel) = df_post_view(i, j, k, MD.c_rev[vel]);
                                }else if(MD.c[vel][1] < 0){
                                    df_view(i, j, k, vel) = df_post_view(i, 16, k, MD.c_rev[vel]);
                                }
                            }
                        }
                    }
                }

                for (int vel = 0; vel < Nvel; vel++)
                {
                    if (norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0)
                    {
                        //pokud miri dovnitr

                        int dx;

                        dx = i - 1;

                        df_view(i, j, k, vel) = cs * df_post_view(dx, j, k, vel) + (1 - cs) *
                            df_post_view(i, j, k, vel);
                    }
                }

                RealType rho = computeRho(i, j, k);

                RealType rho2 = rho_view(i, j, k);

                RealType ux = computeUx(i, j, k, rho2);
                RealType uy = computeUy(i, j, k, rho2);
                RealType uz = computeUz(i, j, k, rho2);

                for (int vel = 0; vel < Nvel; vel++)
                {
                    df_view(i, j, k, vel) = df_view(i, j, k, vel) - f_equilibrium_defined(ux, uy, uz, rho, vel) +
                        f_equilibrium_defined(ux, uy, uz, 1.f, vel);
                }
            }
            /*else if (norm == xm)
            {
                for (int vel = 0; vel < Nvel; vel++)
                {
                    if (norm.x() * MD.c[vel][0] + norm.y() * MD.c[vel][1] + norm.z() * MD.c[vel][2] < 0)
                    {
                        //pokud miri dovnitr

                        int dx;

                        dx = i + 1;

                        df_view(i, j, k, vel) = cs * df_post_view(dx, j, k, vel) + (1 - cs) *
                            df_post_view(i, j, k, vel);
                    }
                }

                RealType rho = computeRho(i, j, k);
                RealType ux = computeUx(i, j, k, density);
                RealType uy = computeUy(i, j, k, density);
                RealType uz = computeUz(i, j, k, density);

                for (int vel = 0; vel < Nvel; vel++)
                {
                    df_view(i, j, k, vel) = df_view(i, j, k, vel) - f_equilibrium_defined(ux, uy, uz, rho, vel) +
                        f_equilibrium_defined(ux, uy, uz, 1.f, vel);
                }
            }*/
            else
            {
                printf("Not yet supported outlet.\n");
            }
        };


        parallelFor<DeviceType>(0, Constants->outlet_num, bb_outlet);
    }
};

#endif
