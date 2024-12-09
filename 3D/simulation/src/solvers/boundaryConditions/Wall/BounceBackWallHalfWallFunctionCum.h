//
// Created by stloufra on 10/30/23.
//

#ifndef BOUNCEBACKWALLHALFWALLFUNC_H
#define BOUNCEBACKWALLHALFWALLFUNC_H

#include "../../../traits/LBMTraits.h"

template <typename MODELDATA>
struct BounceBackWallHalfWallFunc
{
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using TensorType = LBMTraits::TensorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;


    static void slipVelocity(LBMDataPointer& Data, LBMConstantsPointer& Constants)
    {
        auto wall_view = Data->meshBoundaryWall.getView(); //defined in mesher

        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();
        auto mesh_view = Data->meshFluid.getView();

        auto rho_view = Data->rho.getView();

        auto ux_view = Data->ux.getView();
        auto uy_view = Data->uy.getView();
        auto uz_view = Data->uz.getView();

        auto omega_view = Data->omega.getView();

        const auto Nvel = Constants->Nvel;

        MODELDATA MD;

        const RealType nu = Constants->ny; // Kinematic viscosity (adjust for your case)
        const RealType u_inf = Constants->U_inf; // Kinematic viscosity (adjust for your case)
        const RealType h = 0.5f; //TODO: do variable

        const auto cs2 = Constants->cs2;
        const auto Ct = Constants->Ct;

        //Constants LogLaw
        const RealType kappa = 0.41f; //(1.0 / kappa) * std::log(yPlus) + B
        const RealType B = 5.0f;

        //NewtonMethod
        const int maxIter = 100;
        const RealType tol = 1e-5;


        auto NewtonMethod= [=]__cuda_callable__( RealType& tau_w, RealType& u_tau, RealType& ue_magnitude, RealType rho)mutable
        {
            for (int iter = 0; iter < maxIter; ++iter)
            {
                u_tau = std::sqrt(tau_w / rho);
                RealType y_plus = u_tau * h / nu;

                RealType u_plus;
                if (y_plus < 5.f){
                    u_plus = y_plus;
                }
                else{
                    u_plus = (1.f / kappa) * std::log(y_plus) + B;
                }

                RealType residual = u_plus * u_tau - ue_magnitude; //u+=u/u_tau

                if (std::abs(residual) < tol){
                    break;
                }

                //dR/dTau = dR/du_tau du_tau/dtau_w + dR/du+ du+/dy+ dy+/du_tau du_tau/dtau_w
                RealType du_tau_d_tau_w = 0.5f / std::sqrt(tau_w * rho);
                RealType dy_plus_d_tau_w = (1.f / nu) * du_tau_d_tau_w;

                RealType du_plus_d_y_plus = (y_plus < 5.f)
                                                ? 1.f
                                                : (1.f / (kappa * y_plus));

                RealType du_plus_d_tau_w = du_plus_d_y_plus * dy_plus_d_tau_w;

                RealType d_residual_d_tau_w = u_plus * du_tau_d_tau_w + du_plus_d_tau_w * u_tau;

                tau_w -= residual / d_residual_d_tau_w; //Newton x_n+1 = x_n - f(x_n)/f'(x_n)


                if (tau_w < 1e-12) tau_w = 1e-12; // prevent 0 division
            }

        };

        auto stressTensor = [=] __cuda_callable__( int const x, int const y, int const z, RealType const ux, RealType const uy, RealType const uz, RealType const rho, TensorType& tau) mutable
        {

            RealType f_0 = df_view(x,y,z, 0);
            RealType f_1 = df_view(x,y,z, 1);
            RealType f_2 = df_view(x,y,z, 2);
            RealType f_3 = df_view(x,y,z, 3);
            RealType f_4 = df_view(x,y,z, 4);
            RealType f_5 = df_view(x,y,z, 5);
            RealType f_6 = df_view(x,y,z, 6);

            RealType f_7 = df_view(x,y,z, 7);
            RealType f_8 = df_view(x,y,z, 8);
            RealType f_9 = df_view(x,y,z, 9);
            RealType f_10 = df_view(x,y,z, 10);
            RealType f_11 = df_view(x,y,z, 11);
            RealType f_12 = df_view(x,y,z, 12);
            RealType f_13 = df_view(x,y,z, 13);
            RealType f_14 = df_view(x,y,z, 14);
            RealType f_15 = df_view(x,y,z, 15);
            RealType f_16 = df_view(x,y,z, 16);
            RealType f_17 = df_view(x,y,z, 17);
            RealType f_18 = df_view(x,y,z, 18);

            RealType f_19 = df_view(x,y,z, 19);
            RealType f_20 = df_view(x,y,z, 20);
            RealType f_21 = df_view(x,y,z, 21);
            RealType f_22 = df_view(x,y,z, 22);
            RealType f_23 = df_view(x,y,z, 23);
            RealType f_24 = df_view(x,y,z, 24);
            RealType f_25 = df_view(x,y,z, 25);
            RealType f_26 = df_view(x,y,z, 26);

            const RealType omg = omega_view(x,y,z);

            //------------------------------------------------------------------------------------
            //--------------------------- TRANSFORM TO CENTRAL MOMENTS ---------------------------
            //------------------------------------------------------------------------------------

            //Eq Geiger clanek 2015(43)
            //first part of the central moments transformation
            const RealType k_aa0 = (f_21 + f_25) + f_11;
            const RealType k_ab0 = (f_8 + f_10) + f_2;
            const RealType k_ac0 = (f_24 + f_19) + f_15;
            const RealType k_ba0 = (f_14 + f_18) + f_5;
            const RealType k_bb0 = (f_4 + f_3) + f_0;
            const RealType k_bc0 = (f_17 + f_13) + f_6;
            const RealType k_ca0 = (f_20 + f_23) + f_16;
            const RealType k_cb0 = (f_9 + f_7) + f_1;
            const RealType k_cc0 = (f_26 + f_22) + f_12;

            const RealType k_aa1 = (f_21 - f_25) - uz * k_aa0;
            const RealType k_ab1 = (f_8 - f_10) - uz * k_ab0;
            const RealType k_ac1 = (f_24 - f_19) - uz * k_ac0;
            const RealType k_ba1 = (f_14 - f_18) - uz * k_ba0;
            const RealType k_bb1 = (f_4 - f_3) - uz * k_bb0;
            const RealType k_bc1 = (f_17 - f_13) - uz * k_bc0;
            const RealType k_ca1 = (f_20 - f_23) - uz * k_ca0;
            const RealType k_cb1 = (f_9 - f_7) - uz * k_cb0;
            const RealType k_cc1 = (f_26 - f_22) - uz * k_cc0;

            const RealType k_aa2 = (f_21 + f_25) - 2.f * uz * (f_21 - f_25) + uz * uz * k_aa0;
            const RealType k_ab2 = (f_8 + f_10) - 2.f * uz * (f_8 - f_10) + uz * uz * k_ab0;
            const RealType k_ac2 = (f_24 + f_19) - 2.f * uz * (f_24 - f_19) + uz * uz * k_ac0;
            const RealType k_ba2 = (f_14 + f_18) - 2.f * uz * (f_14 - f_18) + uz * uz * k_ba0;
            const RealType k_bb2 = (f_4 + f_3) - 2.f * uz * (f_4 - f_3) + uz * uz * k_bb0;
            const RealType k_bc2 = (f_17 + f_13) - 2.f * uz * (f_17 - f_13) + uz * uz * k_bc0;
            const RealType k_ca2 = (f_20 + f_23) - 2.f * uz * (f_20 - f_23) + uz * uz * k_ca0;
            const RealType k_cb2 = (f_9 + f_7) - 2.f * uz * (f_9 - f_7) + uz * uz * k_cb0;
            const RealType k_cc2 = (f_26 + f_22) - 2.f * uz * (f_26 - f_22) + uz * uz * k_cc0;

            //Eq Geiger clanek 2015(44)
            //second part of the central moments transformation
            const RealType k_a00 = (k_ac0 + k_aa0) + k_ab0;
            const RealType k_b00 = (k_bc0 + k_ba0) + k_bb0;
            const RealType k_c00 = (k_cc0 + k_ca0) + k_cb0;
            const RealType k_a01 = (k_ac1 + k_aa1) + k_ab1;
            const RealType k_b01 = (k_bc1 + k_ba1) + k_bb1;
            const RealType k_c01 = (k_cc1 + k_ca1) + k_cb1;
            const RealType k_a02 = (k_ac2 + k_aa2) + k_ab2;
            const RealType k_b02 = (k_bc2 + k_ba2) + k_bb2;
            const RealType k_c02 = (k_cc2 + k_ca2) + k_cb2;

            const RealType k_a10 = (k_ac0 - k_aa0) - uy * k_a00;
            const RealType k_b10 = (k_bc0 - k_ba0) - uy * k_b00;
            const RealType k_c10 = (k_cc0 - k_ca0) - uy * k_c00;
            const RealType k_a11 = (k_ac1 - k_aa1) - uy * k_a01;
            const RealType k_b11 = (k_bc1 - k_ba1) - uy * k_b01;
            const RealType k_c11 = (k_cc1 - k_ca1) - uy * k_c01;
            const RealType k_a12 = (k_ac2 - k_aa2) - uy * k_a02;
            const RealType k_b12 = (k_bc2 - k_ba2) - uy * k_b02;
            const RealType k_c12 = (k_cc2 - k_ca2) - uy * k_c02;

            const RealType k_a20 = (k_ac0 + k_aa0) - 2.f * uy * (k_ac0 - k_aa0) + uy * uy * k_a00;
            const RealType k_b20 = (k_bc0 + k_ba0) - 2.f * uy * (k_bc0 - k_ba0) + uy * uy * k_b00;
            const RealType k_c20 = (k_cc0 + k_ca0) - 2.f * uy * (k_cc0 - k_ca0) + uy * uy * k_c00;
            const RealType k_a21 = (k_ac1 + k_aa1) - 2.f * uy * (k_ac1 - k_aa1) + uy * uy * k_a01;
            const RealType k_b21 = (k_bc1 + k_ba1) - 2.f * uy * (k_bc1 - k_ba1) + uy * uy * k_b01;
            const RealType k_c21 = (k_cc1 + k_ca1) - 2.f * uy * (k_cc1 - k_ca1) + uy * uy * k_c01;
            const RealType k_a22 = (k_ac2 + k_aa2) - 2.f * uy * (k_ac2 - k_aa2) + uy * uy * k_a02;
            const RealType k_b22 = (k_bc2 + k_ba2) - 2.f * uy * (k_bc2 - k_ba2) + uy * uy * k_b02;
            const RealType k_c22 = (k_cc2 + k_ca2) - 2.f * uy * (k_cc2 - k_ca2) + uy * uy * k_c02;

            //Eq Geiger clanek 2015(45)
            // third part of the central moments transformation
            const RealType k_000 = (k_c00 + k_a00) + k_b00;
            const RealType k_001 = (k_c01 + k_a01) + k_b01;
            const RealType k_002 = (k_c02 + k_a02) + k_b02;
            const RealType k_010 = (k_c10 + k_a10) + k_b10;
            const RealType k_011 = (k_c11 + k_a11) + k_b11;
            const RealType k_020 = (k_c20 + k_a20) + k_b20;

            const RealType k_101 = (k_c01 - k_a01) - ux * k_001;
            const RealType k_110 = (k_c10 - k_a10) - ux * k_010;

            const RealType k_200 = (k_c00 + k_a00) - 2.f * ux * (k_c00 - k_a00) + ux * ux * k_000;


            tau(0)(0) = k_200;
            tau(0)(1) = k_110;
            tau(0)(2) = k_101;
            tau(1)(0) = k_110;
            tau(1)(1) = k_020;
            tau(1)(2) = k_011;
            tau(2)(0) = k_101;
            tau(2)(1) = k_011;
            tau(2)(2) = k_002;

            tau = tau *nu*(-3.f)*omg;
        };



        auto slipVel= [=]__cuda_callable__(const TNL::Containers::StaticArray<1, int>& nod)
        mutable
        {
            Vertex vert = wall_view[nod.x()].vertex;
            bool reg = wall_view[nod.x()].regular;

            if (reg)
            {
                int x = vert.x;
                int y = vert.y;
                int z = vert.z;

                VectorType n = -wall_view[nod.x()].normal; //vnitrni normala

                RealType ux = ux_view(x, y, z);
                RealType uy = uz_view(x, y, z);
                RealType uz = uz_view(x, y, z);

                RealType rho = rho_view(x, y, z);


                VectorType u(ux, uy, uz); //rychlost v prvnim bode
                VectorType e = u - (n, u) * n;
                e = e / TNL::l2Norm(e); //normalizovany smer rychlosti
                VectorType ue = e*(u, e); //tecna rychlost


                RealType ue_magnitude =(u, e);
                RealType tau_w = 0.0001f; //guess

                RealType u_tau(0);

                NewtonMethod(tau_w, u_tau, ue_magnitude, rho);

                RealType mu = rho*nu;

                RealType u_w = ue_magnitude - tau_w/mu/2.f

                /*RealType Cf = 2.f * u_tau * u_tau / u_inf / u_inf;

                TensorType tau;

                stressTensor(x,y,z, ux, uy, uz, rho, tau);

                //uw = -1/Cf (n^T . tau) . e

                //n^T . tau
                VectorType nT_tau(0.f, 0.f, 0.f);

                for (int i = 0; i < 3; ++i){
                    for (int j = 0; j < 3; ++j){
                        nT_tau(i) = nT_tau(i) + n(j) * tau(j)(i);
                    }
                }

                // (n^T . tau) . e
                RealType nT_tau_dot_e = 0.0;
                for (int i = 0; i < 3; ++i){
                    nT_tau_dot_e += nT_tau[i] * e[i];
                }

                // Compute wall velocity u_w
                RealType u_w = -nT_tau_dot_e / Cf;
                */
                e = e*u_w;

                wall_view[nod.x()].slipVelocity = {e(0),e(1),e(2)};*/
            }
        };

        parallelFor<DeviceType>(0, Constants->wall_num, slipVel);
    }


    static void bounceBackWall(LBMDataPointer& Data, LBMConstantsPointer& Constants)
    {
        auto wall_view = Data->meshBoundaryWall.getView(); //defined in mesher
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();
        const auto Nvel = Constants->Nvel;
        auto mesh_view = Data->meshFluid.getView();

        MODELDATA MD;

        auto bb_wall = [=]


        __cuda_callable__(


        const TNL::Containers::StaticArray<1, int>& i
        )
        mutable
        {
            Vertex vert = wall_view[i.x()].vertex;
            bool reg = wall_view[i.x()].regular;

            if (!reg)
            {
                for (int vel = 0; vel < Nvel; vel++)
                {
                    int dx, dy, dz;

                    dx = vert.x - MD.c[vel][0];
                    dy = vert.y - MD.c[vel][1];
                    dz = vert.z - MD.c[vel][2];

                    if (mesh_view(dx, dy, dz) == 0)
                    {
                        df_view(vert.x, vert.y, vert.z, vel) = df_post_view(vert.x, vert.y, vert.z, MD.c_rev[vel]);
                    }
                }
            }
            else
            {
                VectorType e = wall_view[i.x()].slipVelocity;

                for (int vel = 0; vel < Nvel; vel++)
                {
                    int dx, dy, dz;

                    dx = vert.x - MD.c[vel][0];
                    dy = vert.y - MD.c[vel][1];
                    dz = vert.z - MD.c[vel][2];

                    if (mesh_view(dx, dy, dz) == 0)
                    {
                        RealType W = 6 * MD.weight[vel]*(e(0) * MD.c[vel][0] + e(1) * MD.c[vel][1] + e(2) * MD.c[vel][2]);
                        df_view(vert.x, vert.y, vert.z, vel) = df_post_view(vert.x, vert.y, vert.z, MD.c_rev[vel]) - W;
                    }
                }
            }
        };

        parallelFor<DeviceType>(0, Constants->wall_num, bb_wall);
    }
};


/*
 *
 *
 *
 * auto f_equilibrium = [=]__cuda_callable__(
 *                RealType ux,
 *                RealType uy,
 *                RealType uz,
 *                RealType rho,
 *                const int& velo) mutable{
 *           RealType uc, u2;
 *
 *           uc = MD.c[velo][0] * ux
 *               + MD.c[velo][1] * uy
 *               + MD.c[velo][2] * uz;
 *
 *           u2 = ux * ux + uy * uy +
 *               uz * uz;
 *
 *           return MD.weight[velo] * rho * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
 *       };
 *TensorType pi(0.f);
 *
 *for (int vel = 0; vel < Nvel; vel++){
 *RealType feq = f_equilibrium(ux, uy, uz, rho, vel);
 *RealType df = df_view(x, y, z, vel);
 *
 *pi(0)(0) += MD.c[vel][0] * MD.c[vel][0] * (df - feq);
 *pi(0)(1) += MD.c[vel][0] * MD.c[vel][1] * (df - feq);
 *pi(0)(2) += MD.c[vel][0] * MD.c[vel][2] * (df - feq);
 *
 *pi(1)(0) += MD.c[vel][1] * MD.c[vel][0] * (df - feq);
 *pi(1)(1) += MD.c[vel][1] * MD.c[vel][1] * (df - feq);
 *pi(1)(2) += MD.c[vel][1] * MD.c[vel][2] * (df - feq);
 *
 *pi(2)(0) += MD.c[vel][2] * MD.c[vel][0] * (df - feq);
 *pi(2)(1) += MD.c[vel][2] * MD.c[vel][1] * (df - feq);
 *pi(2)(2) += MD.c[vel][2] * MD.c[vel][2] * (df - feq);
 *}
 *
 *RealType StrainRateConst = omega_view(x, y, z) / (-2.f * cs2 * rho * Ct);
 *
 *pi = pi * StrainRateConst; // strain tensor
 *
 *TensorType tau;
 *
 *RealType trace = pi(0)(0) + pi(1)(1) + pi(2)(2);
 *
 *trace = trace / 3.f;
 *
 *RealType mu = nu * rho;
 *
 *for (int i = 0; i < 3; ++i){
 *for (int j = 0; j < 3; ++j){
 *RealType isotropicPart = (i == j) ? trace : 0.0; // Only for diagonal terms
 *tau(i)(j) = 2.0 * mu * (pi(i)(j) - isotropicPart);
 *}
 *}
*/

#endif
