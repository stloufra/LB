//
// Created by stloufra on 10/30/23.
//

#ifndef COLLISIONCUMD3Q27LAMINAR_H
#define COLLISIONCUMD3Q27LAMINAR_H

#include "../../traits/LBMTraits.h"

template<typename MODELDATA>
struct CollisionCumD3Q27Laminar {

    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void collision(LBMDataPointer &Data, LBMConstantsPointer &Constants) {
        auto rho_view = Data->rho.getView();

        auto ux_view = Data->ux.getView();
        auto uy_view = Data->uy.getView();
        auto uz_view = Data->uz.getView();

        auto mesh_view = Data->meshFluid.getView();
        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();


        auto omg = Constants->omega;


        MODELDATA MD;

        auto cum = [=]
        __cuda_callable__(
        const int i,
        const int j,
        const int k)
        mutable
        {
            // M.Geiger 2015

            RealType f_0 = df_view(i, j, k, 0);
            RealType f_1 = df_view(i, j, k, 1);
            RealType f_2 = df_view(i, j, k, 2);
            RealType f_3 = df_view(i, j, k, 3);
            RealType f_4 = df_view(i, j, k, 4);
            RealType f_5 = df_view(i, j, k, 5);
            RealType f_6 = df_view(i, j, k, 6);

            RealType f_7 = df_view(i, j, k, 7);
            RealType f_8 = df_view(i, j, k, 8);
            RealType f_9 = df_view(i, j, k, 9);
            RealType f_10 = df_view(i, j, k, 10);
            RealType f_11 = df_view(i, j, k, 11);
            RealType f_12 = df_view(i, j, k, 12);
            RealType f_13 = df_view(i, j, k, 13);
            RealType f_14 = df_view(i, j, k, 14);
            RealType f_15 = df_view(i, j, k, 15);
            RealType f_16 = df_view(i, j, k, 16);
            RealType f_17 = df_view(i, j, k, 17);
            RealType f_18 = df_view(i, j, k, 18);

            RealType f_19 = df_view(i, j, k, 19);
            RealType f_20 = df_view(i, j, k, 20);
            RealType f_21 = df_view(i, j, k, 21);
            RealType f_22 = df_view(i, j, k, 22);
            RealType f_23 = df_view(i, j, k, 23);
            RealType f_24 = df_view(i, j, k, 24);
            RealType f_25 = df_view(i, j, k, 25);
            RealType f_26 = df_view(i, j, k, 26);

            const RealType ux = ux_view(i, j, k);
            const RealType uy = uy_view(i, j, k);
            const RealType uz = uz_view(i, j, k);

            const RealType rho = rho_view(i, j, k);

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
            const RealType k_012 = (k_c12 + k_a12) + k_b12;
            const RealType k_020 = (k_c20 + k_a20) + k_b20;
            const RealType k_021 = (k_c21 + k_a21) + k_b21;
            const RealType k_022 = (k_c22 + k_a22) + k_b22;

            const RealType k_100 = (k_c00 - k_a00) - ux * k_000;
            const RealType k_101 = (k_c01 - k_a01) - ux * k_001;
            const RealType k_102 = (k_c02 - k_a02) - ux * k_002;
            const RealType k_110 = (k_c10 - k_a10) - ux * k_010;
            const RealType k_111 = (k_c11 - k_a11) - ux * k_011;
            const RealType k_112 = (k_c12 - k_a12) - ux * k_012;
            const RealType k_120 = (k_c20 - k_a20) - ux * k_020;
            const RealType k_121 = (k_c21 - k_a21) - ux * k_021;
            const RealType k_122 = (k_c22 - k_a22) - ux * k_022;

            const RealType k_200 = (k_c00 + k_a00) - 2.f * ux * (k_c00 - k_a00) + ux * ux * k_000;
            const RealType k_201 = (k_c01 + k_a01) - 2.f * ux * (k_c01 - k_a01) + ux * ux * k_001;
            const RealType k_202 = (k_c02 + k_a02) - 2.f * ux * (k_c02 - k_a02) + ux * ux * k_002;
            const RealType k_210 = (k_c10 + k_a10) - 2.f * ux * (k_c10 - k_a10) + ux * ux * k_010;
            const RealType k_211 = (k_c11 + k_a11) - 2.f * ux * (k_c11 - k_a11) + ux * ux * k_011;
            const RealType k_212 = (k_c12 + k_a12) - 2.f * ux * (k_c12 - k_a12) + ux * ux * k_012;
            const RealType k_220 = (k_c20 + k_a20) - 2.f * ux * (k_c20 - k_a20) + ux * ux * k_020;
            const RealType k_221 = (k_c21 + k_a21) - 2.f * ux * (k_c21 - k_a21) + ux * ux * k_021;
            const RealType k_222 = (k_c22 + k_a22) - 2.f * ux * (k_c22 - k_a22) + ux * ux * k_022;

            //------------------------------------------------------------------------------------
            //------------------------------ CENTRAL MOM. TO CUMULANTS ---------------------------
            //------------------------------------------------------------------------------------


            //-------------Until order 3 proportional to the central moments---------------

            //Eq Geiger clanek 2015(47)
            const RealType C_110 = k_110;
            const RealType C_101 = k_101;
            const RealType C_011 = k_011;

            //Eq Geiger clanek 2015(48)
            const RealType C_200 = k_200;
            const RealType C_020 = k_020;
            const RealType C_002 = k_002;

            //Eq Geiger clanek 2015(49)
            const RealType C_120 = k_120;
            const RealType C_012 = k_012;
            const RealType C_201 = k_201;

            const RealType C_210 = k_210;
            const RealType C_021 = k_021;
            const RealType C_102 = k_102;

            //Eq Geiger clanek 2015(50)
            const RealType C_111 = k_111;

            //Eq Geiger clanek 2015(51)
            const RealType C_211 = k_211 - (k_200 * k_011 + 2.f * k_101 * k_110) / rho;
            const RealType C_121 = k_121 - (k_020 * k_101 + 2.f * k_110 * k_011) / rho;
            const RealType C_112 = k_112 - (k_002 * k_110 + 2.f * k_011 * k_101) / rho;

            //Eq Geiger clanek 2015(52)
            const RealType C_220 = k_220 - (k_020 * k_200 + 2.f * k_110 * k_110) / rho;
            const RealType C_022 = k_022 - (k_002 * k_020 + 2.f * k_011 * k_011) / rho;
            const RealType C_202 = k_202 - (k_200 * k_002 + 2.f * k_101 * k_101) / rho;

            //Eq Geiger clanek 2015(53)
            const RealType C_122 = k_122 - (k_020 * k_102 + k_002 * k_120 + 4.f * k_011 * k_111 +
                2.f * (k_110 * k_012 + k_101 * k_021)) / rho;
            const RealType C_212 = k_212 - (k_002 * k_210 + k_200 * k_012 + 4.f * k_101 * k_111 +
                2.f * (k_011 * k_201 + k_110 * k_102)) / rho;
            const RealType C_221 = k_221 - (k_200 * k_021 + k_020 * k_201 + 4.f * k_110 * k_111 +
                2.f * (k_101 * k_120 + k_011 * k_210)) / rho;
            //Eq Geiger clanek 2015(54)
            const RealType C_222 = k_222 - (4.f * k_111 * k_111 + k_200 * k_022 + k_020 * k_202 + k_002 * k_220 +
                    4.f * (k_011 * k_211 + k_101 * k_121 + k_110 * k_112) +
                    2.f * (k_120 * k_102 + k_210 * k_012 + k_201 * k_021)) / rho +
                (16.f * k_110 * k_101 * k_011 +
                    4.f * (k_101 * k_101 * k_020 + k_011 * k_011 * k_200 + k_110 * k_110 * k_002) +
                    2.f * k_200 * k_020 * k_002) / rho / rho;

            //------------------------------------------------------------------------------------
            // -------------------------------------COLLISION-------------------------------------
            //------------------------------------------------------------------------------------

            //  RELAX RATE Geiger clanek 2015(103) //2017 diff

            const RealType omega1 = omg;
            const RealType omega2 = 1.f;
            const RealType omega3 = 1.f;
            const RealType omega4 = 1.f;
            const RealType omega5 = 1.f;
            const RealType omega6 = 1.f;
            const RealType omega7 = 1.f;
            const RealType omega8 = 1.f;
            const RealType omega9 = 1.f;
            const RealType omega10 = 1.f;

            //Eq Geiger clanek 2015(58)
            const RealType Dxu = -omega1 * 0.5f / rho * (2.f * C_200 - C_020 - C_002) -
                omega2 * 0.5f / rho * (C_200 + C_020 + C_002 - k_000); // -(-1-rho));

            //Eq Geiger clanek 2015(59)
            const RealType Dyv = Dxu + 3.f * omega1 * 0.5f / rho * (C_200 - C_020);

            //Eq Geiger clanek 2015(60)
            const RealType Dzw = Dxu + 3.f * omega1 * 0.5f / rho * (C_200 - C_002);

            //------------------------------------------------------------------------------------

            //Eq Geiger clanek 2015(55)
            const RealType Cs_110 = (1.f - omega1) * C_110;
            //Eq Geiger clanek 2015(56)
            const RealType Cs_101 = (1.f - omega1) * C_101;
            //Eq Geiger clanek 2015(57)
            const RealType Cs_011 = (1.f - omega1) * C_011;

            //---------------------------------------------------------------------------------

            //Eq Geiger clanek 2015(61, 62, 63)
            const RealType Eq61RHS = (1.f - omega1) * (C_200 - C_020) -
                3.f * rho * (1.f - omega1 * 0.5f) * (ux * ux * Dxu - uy * uy * Dyv);
            const RealType Eq64RHS = (1.f - omega1) * (C_200 - C_002) -
                3.f * rho * (1.f - omega1 * 0.5f) * (ux * ux * Dxu - uz * uz * Dzw);
            const RealType Eq65RHS = k_000 * omega2 + (1.f - omega2) * (C_200 + C_020 + C_002) -
                3.f * rho * (1.f - omega2 / 2.f) * (ux * ux * Dxu + uy * uy * Dyv + uz * uz * Dzw);

            const RealType Cs_200 = 1.f / 3.f * (Eq61RHS + Eq64RHS + Eq65RHS);
            const RealType Cs_020 = 1.f / 3.f * (Eq64RHS - 2.f * Eq61RHS + Eq65RHS);
            const RealType Cs_002 = 1.f / 3.f * (Eq61RHS - 2.f * Eq64RHS + Eq65RHS);

            //Eq Geiger clanek 2015(64, 67)
            const RealType Cs_120 = (-C_102 - C_120) * omega3 * 0.5f + (C_102 - C_120) * omega4 * 0.5f + C_120;
            const RealType Cs_102 = (-C_102 - C_120) * omega3 * 0.5f + (-C_102 + C_120) * omega4 * 0.5f + C_102;

            //Eq Geiger clanek 2015(65, 68)
            const RealType Cs_210 = (-C_012 - C_210) * omega3 * 0.5f + (C_012 - C_210) * omega4 * 0.5f + C_210;
            const RealType Cs_012 = (-C_012 - C_210) * omega3 * 0.5f + (-C_012 + C_210) * omega4 * 0.5f + C_012;

            //Eq Geiger clanek 2015(66, 69)
            const RealType Cs_021 = (-C_021 - C_201) * omega3 * 0.5f + (-C_021 + C_201) * omega4 * 0.5f + C_021;
            const RealType Cs_201 = (-C_021 - C_201) * omega3 * 0.5f + (C_021 - C_201) * omega4 * 0.5f + C_201;

            //Eq Geiger clanek 2015(70)
            const RealType Cs_111 = (1.f - omega5) * C_111;

            //Eq Geiger clanek 2015(71, 72, 73) //2017 diff
            const RealType Eq71RHS = (1.f - omega6) * (C_220 - 2.f * C_202 + C_022);
            const RealType Eq72RHS = (1.f - omega6) * (C_220 + C_202 - 2.f * C_022);
            const RealType Eq73RHS = (1.f - omega7) * (C_220 + C_202 + C_022);

            const RealType Cs_220 = 1.f / 3.f * (Eq71RHS + Eq72RHS + Eq73RHS);
            const RealType Cs_202 = 1.f / 3.f * (-Eq71RHS + Eq73RHS);
            const RealType Cs_022 = 1.f / 3.f * (-Eq72RHS + Eq73RHS);

            // Eq 46-48 //2017 diff
            const RealType Cs_211 = (1.f - omega8) * C_211;
            const RealType Cs_121 = (1.f - omega8) * C_121;
            const RealType Cs_112 = (1.f - omega8) * C_112;
            // Eqs 49-52
            const RealType Cs_221 = (1.f - omega9) * C_221;
            const RealType Cs_212 = (1.f - omega9) * C_212;
            const RealType Cs_122 = (1.f - omega9) * C_122;
            const RealType Cs_222 = (1.f - omega10) * C_222;

            //------------------------------------------------------------------------------------
            //------------------------------ CUMULANTS TO CENTRAL MOM. ---------------------------
            //------------------------------------------------------------------------------------

            const RealType ks_000 = k_000;

            // Permutation again

            //Eq Geiger clanek 2015(47) backwards
            const RealType ks_110 = Cs_110;
            const RealType ks_101 = Cs_101;
            const RealType ks_011 = Cs_011;

            //Eq Geiger clanek 2015(48) backwards
            const RealType ks_200 = Cs_200;
            const RealType ks_020 = Cs_020;
            const RealType ks_002 = Cs_002;

            //Eq Geiger clanek 2015(49) backwards
            const RealType ks_120 = Cs_120;
            const RealType ks_012 = Cs_012;
            const RealType ks_201 = Cs_201;

            const RealType ks_210 = Cs_210;
            const RealType ks_021 = Cs_021;
            const RealType ks_102 = Cs_102;

            //Eq Geiger clanek 2015(50) backwards
            const RealType ks_111 = Cs_111;

            //Eq. Geiger clanek 2015(85, 86, 87)
            const RealType ks_100 = -k_100;
            const RealType ks_010 = -k_010;
            const RealType ks_001 = -k_001;

            //Eq. Geiger clanek 2015(81)
            const RealType ks_211 = Cs_211 + (ks_200 * ks_011 + 2.f * ks_101 * ks_110) / rho;
            const RealType ks_121 = Cs_121 + (ks_020 * ks_101 + 2.f * ks_110 * ks_011) / rho;
            const RealType ks_112 = Cs_112 + (ks_002 * ks_110 + 2.f * ks_011 * ks_101) / rho;

            //Eq. Geiger clanek 2015(82)
            const RealType ks_220 = Cs_220 + (ks_020 * ks_200 + 2.f * ks_110 * ks_110) / rho;
            const RealType ks_022 = Cs_022 + (ks_002 * ks_020 + 2.f * ks_011 * ks_011) / rho;
            const RealType ks_202 = Cs_202 + (ks_200 * ks_002 + 2.f * ks_101 * ks_101) / rho;

            //Eq. Geiger clanek 2015(83)
            const RealType ks_122 = Cs_122 + (ks_020 * ks_102 + ks_002 * ks_120 + 4.f * ks_011 * ks_111 +
                2.f * (ks_110 * ks_012 + ks_101 * ks_021)) / rho;
            const RealType ks_212 = Cs_212 + (ks_002 * ks_210 + ks_200 * ks_012 + 4.f * ks_101 * ks_111 +
                2.f * (ks_011 * ks_201 + ks_110 * ks_102)) / rho;
            const RealType ks_221 = Cs_221 + (ks_200 * ks_021 + ks_020 * ks_201 + 4.f * ks_110 * ks_111 +
                2.f * (ks_101 * ks_120 + ks_011 * ks_210)) / rho;
            // gen2nowell.php END

            // Eq. Geiger clanek 2015(84)
            const RealType ks_222 = Cs_222 +
                (4.f * ks_111 * ks_111 + ks_200 * ks_022 + ks_020 * ks_202 + ks_002 * ks_220 +
                    4.f * (ks_011 * ks_211 + ks_101 * ks_121 + ks_110 * ks_112)
                    + 2.f * (ks_120 * ks_102 + ks_210 * ks_012 + ks_201 * ks_021)) / rho
                - (16.f * ks_110 * ks_101 * ks_011 + 4.f * (ks_101 * ks_101 * ks_020 +
                        ks_011 * ks_011 * ks_200 +
                        ks_110 * ks_110 * ks_002) +
                    2.f * ks_200 * ks_020 * ks_002) / rho / rho;



            //------------------------------------------------------------------------------------
            //----------------------- TRANSFORM TO DISTRIBUTION FUNCTION -------------------------
            //------------------------------------------------------------------------------------

            //Eq Geiger clanek 2015(88)
            const RealType ks_b00 = ks_000 * (1.f - ux * ux) - 2.f * ux * ks_100 - ks_200;
            const RealType ks_b01 = ks_001 * (1.f - ux * ux) - 2.f * ux * ks_101 - ks_201;
            const RealType ks_b02 = ks_002 * (1.f - ux * ux) - 2.f * ux * ks_102 - ks_202;
            const RealType ks_b10 = ks_010 * (1.f - ux * ux) - 2.f * ux * ks_110 - ks_210;
            const RealType ks_b11 = ks_011 * (1.f - ux * ux) - 2.f * ux * ks_111 - ks_211;
            const RealType ks_b12 = ks_012 * (1.f - ux * ux) - 2.f * ux * ks_112 - ks_212;
            const RealType ks_b20 = ks_020 * (1.f - ux * ux) - 2.f * ux * ks_120 - ks_220;
            const RealType ks_b21 = ks_021 * (1.f - ux * ux) - 2.f * ux * ks_121 - ks_221;
            const RealType ks_b22 = ks_022 * (1.f - ux * ux) - 2.f * ux * ks_122 - ks_222;

            //Eq  Geiger clanek 2015(89)
            const RealType ks_a00 = (ks_000 * (ux * ux - ux) + ks_100 * (2.f * ux - 1.f) + ks_200) * 0.5f;
            const RealType ks_a01 = (ks_001 * (ux * ux - ux) + ks_101 * (2.f * ux - 1.f) + ks_201) * 0.5f;
            const RealType ks_a02 = (ks_002 * (ux * ux - ux) + ks_102 * (2.f * ux - 1.f) + ks_202) * 0.5f;
            const RealType ks_a10 = (ks_010 * (ux * ux - ux) + ks_110 * (2.f * ux - 1.f) + ks_210) * 0.5f;
            const RealType ks_a11 = (ks_011 * (ux * ux - ux) + ks_111 * (2.f * ux - 1.f) + ks_211) * 0.5f;
            const RealType ks_a12 = (ks_012 * (ux * ux - ux) + ks_112 * (2.f * ux - 1.f) + ks_212) * 0.5f;
            const RealType ks_a20 = (ks_020 * (ux * ux - ux) + ks_120 * (2.f * ux - 1.f) + ks_220) * 0.5f;
            const RealType ks_a21 = (ks_021 * (ux * ux - ux) + ks_121 * (2.f * ux - 1.f) + ks_221) * 0.5f;
            const RealType ks_a22 = (ks_022 * (ux * ux - ux) + ks_122 * (2.f * ux - 1.f) + ks_222) * 0.5f;

            //Eq  Geiger clanek 2015(90)
            const RealType ks_c00 = (ks_000 * (ux * ux + ux) + ks_100 * (2.f * ux + 1.f) + ks_200) * 0.5f;
            const RealType ks_c01 = (ks_001 * (ux * ux + ux) + ks_101 * (2.f * ux + 1.f) + ks_201) * 0.5f;
            const RealType ks_c02 = (ks_002 * (ux * ux + ux) + ks_102 * (2.f * ux + 1.f) + ks_202) * 0.5f;
            const RealType ks_c10 = (ks_010 * (ux * ux + ux) + ks_110 * (2.f * ux + 1.f) + ks_210) * 0.5f;
            const RealType ks_c11 = (ks_011 * (ux * ux + ux) + ks_111 * (2.f * ux + 1.f) + ks_211) * 0.5f;
            const RealType ks_c12 = (ks_012 * (ux * ux + ux) + ks_112 * (2.f * ux + 1.f) + ks_212) * 0.5f;
            const RealType ks_c20 = (ks_020 * (ux * ux + ux) + ks_120 * (2.f * ux + 1.f) + ks_220) * 0.5f;
            const RealType ks_c21 = (ks_021 * (ux * ux + ux) + ks_121 * (2.f * ux + 1.f) + ks_221) * 0.5f;
            const RealType ks_c22 = (ks_022 * (ux * ux + ux) + ks_122 * (2.f * ux + 1.f) + ks_222) * 0.5f;

            //Eq Geiger clanek 2015(91)
            const RealType ks_ab0 = ks_a00 * (1.f - uy * uy) - 2.f * uy * ks_a10 - ks_a20;
            const RealType ks_ab1 = ks_a01 * (1.f - uy * uy) - 2.f * uy * ks_a11 - ks_a21;
            const RealType ks_ab2 = ks_a02 * (1.f - uy * uy) - 2.f * uy * ks_a12 - ks_a22;
            const RealType ks_bb0 = ks_b00 * (1.f - uy * uy) - 2.f * uy * ks_b10 - ks_b20;
            const RealType ks_bb1 = ks_b01 * (1.f - uy * uy) - 2.f * uy * ks_b11 - ks_b21;
            const RealType ks_bb2 = ks_b02 * (1.f - uy * uy) - 2.f * uy * ks_b12 - ks_b22;
            const RealType ks_cb0 = ks_c00 * (1.f - uy * uy) - 2.f * uy * ks_c10 - ks_c20;
            const RealType ks_cb1 = ks_c01 * (1.f - uy * uy) - 2.f * uy * ks_c11 - ks_c21;
            const RealType ks_cb2 = ks_c02 * (1.f - uy * uy) - 2.f * uy * ks_c12 - ks_c22;

            //Eq  Geiger clanek 2015(92)
            const RealType ks_aa0 = (ks_a00 * (uy * uy - uy) + ks_a10 * (2.f * uy - 1.f) + ks_a20) * 0.5f;
            const RealType ks_aa1 = (ks_a01 * (uy * uy - uy) + ks_a11 * (2.f * uy - 1.f) + ks_a21) * 0.5f;
            const RealType ks_aa2 = (ks_a02 * (uy * uy - uy) + ks_a12 * (2.f * uy - 1.f) + ks_a22) * 0.5f;
            const RealType ks_ba0 = (ks_b00 * (uy * uy - uy) + ks_b10 * (2.f * uy - 1.f) + ks_b20) * 0.5f;
            const RealType ks_ba1 = (ks_b01 * (uy * uy - uy) + ks_b11 * (2.f * uy - 1.f) + ks_b21) * 0.5f;
            const RealType ks_ba2 = (ks_b02 * (uy * uy - uy) + ks_b12 * (2.f * uy - 1.f) + ks_b22) * 0.5f;
            const RealType ks_ca0 = (ks_c00 * (uy * uy - uy) + ks_c10 * (2.f * uy - 1.f) + ks_c20) * 0.5f;
            const RealType ks_ca1 = (ks_c01 * (uy * uy - uy) + ks_c11 * (2.f * uy - 1.f) + ks_c21) * 0.5f;
            const RealType ks_ca2 = (ks_c02 * (uy * uy - uy) + ks_c12 * (2.f * uy - 1.f) + ks_c22) * 0.5f;

            //Eq Geiger clanek 2015(93)
            const RealType ks_ac0 = (ks_a00 * (uy * uy + uy) + ks_a10 * (2.f * uy + 1.f) + ks_a20) * 0.5f;
            const RealType ks_ac1 = (ks_a01 * (uy * uy + uy) + ks_a11 * (2.f * uy + 1.f) + ks_a21) * 0.5f;
            const RealType ks_ac2 = (ks_a02 * (uy * uy + uy) + ks_a12 * (2.f * uy + 1.f) + ks_a22) * 0.5f;
            const RealType ks_bc0 = (ks_b00 * (uy * uy + uy) + ks_b10 * (2.f * uy + 1.f) + ks_b20) * 0.5f;
            const RealType ks_bc1 = (ks_b01 * (uy * uy + uy) + ks_b11 * (2.f * uy + 1.f) + ks_b21) * 0.5f;
            const RealType ks_bc2 = (ks_b02 * (uy * uy + uy) + ks_b12 * (2.f * uy + 1.f) + ks_b22) * 0.5f;
            const RealType ks_cc0 = (ks_c00 * (uy * uy + uy) + ks_c10 * (2.f * uy + 1.f) + ks_c20) * 0.5f;
            const RealType ks_cc1 = (ks_c01 * (uy * uy + uy) + ks_c11 * (2.f * uy + 1.f) + ks_c21) * 0.5f;
            const RealType ks_cc2 = (ks_c02 * (uy * uy + uy) + ks_c12 * (2.f * uy + 1.f) + ks_c22) * 0.5f;

            //Eq Geiger clanek 2015(94)
            f_11 = ks_aa0 * (1.f - uz * uz) - 2.f * uz * ks_aa1 - ks_aa2;
            f_2 = ks_ab0 * (1.f - uz * uz) - 2.f * uz * ks_ab1 - ks_ab2;
            f_15 = ks_ac0 * (1.f - uz * uz) - 2.f * uz * ks_ac1 - ks_ac2;
            f_5 = ks_ba0 * (1.f - uz * uz) - 2.f * uz * ks_ba1 - ks_ba2;
            f_0 = ks_bb0 * (1.f - uz * uz) - 2.f * uz * ks_bb1 - ks_bb2;
            f_6 = ks_bc0 * (1.f - uz * uz) - 2.f * uz * ks_bc1 - ks_bc2;
            f_16 = ks_ca0 * (1.f - uz * uz) - 2.f * uz * ks_ca1 - ks_ca2;
            f_1 = ks_cb0 * (1.f - uz * uz) - 2.f * uz * ks_cb1 - ks_cb2;
            f_12 = ks_cc0 * (1.f - uz * uz) - 2.f * uz * ks_cc1 - ks_cc2;

            //Eq  Geiger clanek 2015(95)
            f_25 = (ks_aa0 * (uz * uz - uz) + ks_aa1 * (2.f * uz - 1.f) + ks_aa2) * 0.5f;
            f_10 = (ks_ab0 * (uz * uz - uz) + ks_ab1 * (2.f * uz - 1.f) + ks_ab2) * 0.5f;
            f_19 = (ks_ac0 * (uz * uz - uz) + ks_ac1 * (2.f * uz - 1.f) + ks_ac2) * 0.5f;
            f_18 = (ks_ba0 * (uz * uz - uz) + ks_ba1 * (2.f * uz - 1.f) + ks_ba2) * 0.5f;
            f_3 = (ks_bb0 * (uz * uz - uz) + ks_bb1 * (2.f * uz - 1.f) + ks_bb2) * 0.5f;
            f_13 = (ks_bc0 * (uz * uz - uz) + ks_bc1 * (2.f * uz - 1.f) + ks_bc2) * 0.5f;
            f_23 = (ks_ca0 * (uz * uz - uz) + ks_ca1 * (2.f * uz - 1.f) + ks_ca2) * 0.5f;
            f_7 = (ks_cb0 * (uz * uz - uz) + ks_cb1 * (2.f * uz - 1.f) + ks_cb2) * 0.5f;
            f_22 = (ks_cc0 * (uz * uz - uz) + ks_cc1 * (2.f * uz - 1.f) + ks_cc2) * 0.5f;

            //Eq  Geiger clanek 2015(96)
            f_21 = (ks_aa0 * (uz * uz + uz) + ks_aa1 * (2.f * uz + 1.f) + ks_aa2) * 0.5f;
            f_8 = (ks_ab0 * (uz * uz + uz) + ks_ab1 * (2.f * uz + 1.f) + ks_ab2) * 0.5f;
            f_24 = (ks_ac0 * (uz * uz + uz) + ks_ac1 * (2.f * uz + 1.f) + ks_ac2) * 0.5f;
            f_14 = (ks_ba0 * (uz * uz + uz) + ks_ba1 * (2.f * uz + 1.f) + ks_ba2) * 0.5f;
            f_4 = (ks_bb0 * (uz * uz + uz) + ks_bb1 * (2.f * uz + 1.f) + ks_bb2) * 0.5f;
            f_17 = (ks_bc0 * (uz * uz + uz) + ks_bc1 * (2.f * uz + 1.f) + ks_bc2) * 0.5f;
            f_20 = (ks_ca0 * (uz * uz + uz) + ks_ca1 * (2.f * uz + 1.f) + ks_ca2) * 0.5f;
            f_9 = (ks_cb0 * (uz * uz + uz) + ks_cb1 * (2.f * uz + 1.f) + ks_cb2) * 0.5f;
            f_26 = (ks_cc0 * (uz * uz + uz) + ks_cc1 * (2.f * uz + 1.f) + ks_cc2) * 0.5f;

            // back to df_post
            df_post_view(i, j, k, 0) = f_0; // 0
            df_post_view(i, j, k, 1) = f_1; // 1
            df_post_view(i, j, k, 2) = f_2; // 2
            df_post_view(i, j, k, 6) = f_6; // 6
            df_post_view(i, j, k, 5) = f_5; // 5
            df_post_view(i, j, k, 4) = f_4; // 4
            df_post_view(i, j, k, 3) = f_3; // 3

            df_post_view(i, j, k, 12) = f_12; // 12
            df_post_view(i, j, k, 11) = f_11; // 11
            df_post_view(i, j, k, 16) = f_16; // 16
            df_post_view(i, j, k, 15) = f_15; // 15
            df_post_view(i, j, k, 9) = f_9; // 9
            df_post_view(i, j, k, 10) = f_10; // 10
            df_post_view(i, j, k, 7) = f_7; // 7
            df_post_view(i, j, k, 8) = f_8; // 8
            df_post_view(i, j, k, 17) = f_17; // 17
            df_post_view(i, j, k, 18) = f_18; // 18
            df_post_view(i, j, k, 13) = f_13; // 13
            df_post_view(i, j, k, 14) = f_14; // 14

            df_post_view(i, j, k, 26) = f_26; // 26
            df_post_view(i, j, k, 25) = f_25; // 25
            df_post_view(i, j, k, 22) = f_22; // 22
            df_post_view(i, j, k, 21) = f_21; // 21
            df_post_view(i, j, k, 20) = f_20; // 20
            df_post_view(i, j, k, 19) = f_19; // 19
            df_post_view(i, j, k, 23) = f_23; // 23
            df_post_view(i, j, k, 24) = f_24; // 24
        };


        auto coll = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i )mutable
        {

            if (mesh_view(i.x(), i.y(), i.z()) != 0) {

                cum(i.x(), i.y(), i.z());
            }
        };


        TNL::Containers::StaticArray < 3, int > begin{1, 1, 1};
        TNL::Containers::StaticArray < 3, int > end{Constants->dimX_int-1, Constants->dimY_int-1, Constants->dimZ_int-1};
        parallelFor<DeviceType>(begin, end, coll);
    }


};


#endif
