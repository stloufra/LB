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

        auto f_equilibrium = [=]__cuda_callable__(
                 RealType ux,
                 RealType uy,
                 RealType uz,
                 RealType rho,
                 const int& velo)
        mutable
        {
            RealType uc, u2;

            uc = MD.c[velo][0] * ux
                + MD.c[velo][1] * uy
                + MD.c[velo][2] * uz;

            u2 = ux * ux + uy * uy +
                uz * uz;


            return MD.weight[velo] * rho * (1.f + 3.f * uc + 4.5f * uc * uc - 1.5f * u2);
        };

        auto bb_wall = [=]__cuda_callable__(const TNL::Containers::StaticArray<1, int>& nod)
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
                VectorType e = u - (u, n) * n;
                e = e / TNL::l2Norm(e); //normalizovany smer rychlosti
                VectorType ue = (u, e); //tecna rychlost


                RealType ue_magnitude =TNL::l2Norm(ue);
                RealType tau_w = 0.0001f; //guess

                RealType u_tau;

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

                RealType Cf = 2.f * u_tau * u_tau / u_inf / u_inf;
                TensorType pi(0.f);

                for (int vel = 0; vel < Nvel; vel++){
                    RealType feq = f_equilibrium(ux, uy, uz, rho, vel);
                    RealType df = df_view(x, y, z, vel);

                    pi(0)(0) += MD.c[vel][0] * MD.c[vel][0] * (df - feq);
                    pi(0)(1) += MD.c[vel][0] * MD.c[vel][1] * (df - feq);
                    pi(0)(2) += MD.c[vel][0] * MD.c[vel][2] * (df - feq);

                    pi(1)(0) += MD.c[vel][1] * MD.c[vel][0] * (df - feq);
                    pi(1)(1) += MD.c[vel][1] * MD.c[vel][1] * (df - feq);
                    pi(1)(2) += MD.c[vel][1] * MD.c[vel][2] * (df - feq);

                    pi(2)(0) += MD.c[vel][2] * MD.c[vel][0] * (df - feq);
                    pi(2)(1) += MD.c[vel][2] * MD.c[vel][1] * (df - feq);
                    pi(2)(2) += MD.c[vel][2] * MD.c[vel][2] * (df - feq);
                }

                RealType StrainRateConst = omega_view(x, y, z) / (-2 * cs2 * rho * Ct);

                pi = pi * StrainRateConst; // strain tensor

                TensorType tau;

                RealType trace = pi(0)(0) + pi(1)(1) + pi(2)(2);

                trace = trace / 3;

                RealType mu = nu * rho;

                for (int i = 0; i < 3; ++i){
                    for (int j = 0; j < 3; ++j){
                        RealType isotropicPart = (i == j) ? trace : 0.0; // Only for diagonal terms
                        tau(i)(j) = 2.0 * mu * (pi(i)(j) - isotropicPart);
                    }
                }

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

                e = e*u_w;

                wall_view[nod.x()].slipVelocity = e;
            }
        };

        parallelFor<DeviceType>(0, Constants->wall_num, bb_wall);
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
                        RealType W = 6 * MD.weight[
                            vel](e(0) * MD.c[vel][0] + e(1) * MD.c[vel][1] + e(2) * MD.c[vel][2]);
                        df_view(vert.x, vert.y, vert.z, vel) = df_post_view(vert.x, vert.y, vert.z, MD.c_rev[vel]) - W;
                    }
                }
            }
        };

        parallelFor<DeviceType>(0, Constants->wall_num, bb_wall);
    }
};

#endif
