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
        auto u_view = Data->u.getView();
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

            // for easier implementation (TNL-LBM naming convention with their dynamics)
            // implementation of M.Geiger 2015

            RealType f_zzz = df_view(i, j, k, 0);
            RealType f_pzz = df_view(i, j, k, 1);
            RealType f_mzz = df_view(i, j, k, 2);
            RealType f_zpz = df_view(i, j, k, 6);
            RealType f_zmz = df_view(i, j, k, 5);
            RealType f_zzp = df_view(i, j, k, 4);
            RealType f_zzm = df_view(i, j, k, 3);

            RealType f_ppz = df_view(i, j, k, 12);
            RealType f_mmz = df_view(i, j, k, 11);
            RealType f_pmz = df_view(i, j, k, 16);
            RealType f_mpz = df_view(i, j, k, 15);
            RealType f_pzp = df_view(i, j, k, 9);
            RealType f_mzm = df_view(i, j, k, 10);
            RealType f_pzm = df_view(i, j, k, 7);
            RealType f_mzp = df_view(i, j, k, 8);
            RealType f_zpp = df_view(i, j, k, 17);
            RealType f_zmm = df_view(i, j, k, 18);
            RealType f_zpm = df_view(i, j, k, 13);
            RealType f_zmp = df_view(i, j, k, 14);

            RealType f_ppp = df_view(i, j, k, 26);
            RealType f_mmm = df_view(i, j, k, 25);
            RealType f_ppm = df_view(i, j, k, 22);
            RealType f_mmp = df_view(i, j, k, 21);
            RealType f_pmp = df_view(i, j, k, 20);
            RealType f_mpm = df_view(i, j, k, 19);
            RealType f_pmm = df_view(i, j, k, 23);
            RealType f_mpp = df_view(i, j, k, 24);

            const RealType ux = u_view(i, j, k).x();
            const RealType uy = u_view(i, j, k).y();
            const RealType uz = u_view(i, j, k).z();

            const RealType rho = rho_view(i, j, k);

            //------------------------------------------------------------------------------------
            //--------------------------- TRANSFORM TO CENTRAL MOMENTS ---------------------------
            //------------------------------------------------------------------------------------

            //Eq G2015(43)
            //first part of the central moments transformation
            const RealType k_mm0 = (f_mmp + f_mmm) + f_mmz;
            const RealType k_mz0 = (f_mzp + f_mzm) + f_mzz;
            const RealType k_mp0 = (f_mpp + f_mpm) + f_mpz;
            const RealType k_zm0 = (f_zmp + f_zmm) + f_zmz;
            const RealType k_zz0 = (f_zzp + f_zzm) + f_zzz;
            const RealType k_zp0 = (f_zpp + f_zpm) + f_zpz;
            const RealType k_pm0 = (f_pmp + f_pmm) + f_pmz;
            const RealType k_pz0 = (f_pzp + f_pzm) + f_pzz;
            const RealType k_pp0 = (f_ppp + f_ppm) + f_ppz;

            const RealType k_mm1 = (f_mmp - f_mmm) - uz * k_mm0;
            const RealType k_mz1 = (f_mzp - f_mzm) - uz * k_mz0;
            const RealType k_mp1 = (f_mpp - f_mpm) - uz * k_mp0;
            const RealType k_zm1 = (f_zmp - f_zmm) - uz * k_zm0;
            const RealType k_zz1 = (f_zzp - f_zzm) - uz * k_zz0;
            const RealType k_zp1 = (f_zpp - f_zpm) - uz * k_zp0;
            const RealType k_pm1 = (f_pmp - f_pmm) - uz * k_pm0;
            const RealType k_pz1 = (f_pzp - f_pzm) - uz * k_pz0;
            const RealType k_pp1 = (f_ppp - f_ppm) - uz * k_pp0;

            const RealType k_mm2 = (f_mmp + f_mmm) - 2.f * uz * (f_mmp - f_mmm) + uz * uz * k_mm0;
            const RealType k_mz2 = (f_mzp + f_mzm) - 2.f * uz * (f_mzp - f_mzm) + uz * uz * k_mz0;
            const RealType k_mp2 = (f_mpp + f_mpm) - 2.f * uz * (f_mpp - f_mpm) + uz * uz * k_mp0;
            const RealType k_zm2 = (f_zmp + f_zmm) - 2.f * uz * (f_zmp - f_zmm) + uz * uz * k_zm0;
            const RealType k_zz2 = (f_zzp + f_zzm) - 2.f * uz * (f_zzp - f_zzm) + uz * uz * k_zz0;
            const RealType k_zp2 = (f_zpp + f_zpm) - 2.f * uz * (f_zpp - f_zpm) + uz * uz * k_zp0;
            const RealType k_pm2 = (f_pmp + f_pmm) - 2.f * uz * (f_pmp - f_pmm) + uz * uz * k_pm0;
            const RealType k_pz2 = (f_pzp + f_pzm) - 2.f * uz * (f_pzp - f_pzm) + uz * uz * k_pz0;
            const RealType k_pp2 = (f_ppp + f_ppm) - 2.f * uz * (f_ppp - f_ppm) + uz * uz * k_pp0;

            //Eq G2015(44)
            //second part of the central moments transformation
            const RealType k_m00 = (k_mp0 + k_mm0) + k_mz0;
            const RealType k_z00 = (k_zp0 + k_zm0) + k_zz0;
            const RealType k_p00 = (k_pp0 + k_pm0) + k_pz0;
            const RealType k_m01 = (k_mp1 + k_mm1) + k_mz1;
            const RealType k_z01 = (k_zp1 + k_zm1) + k_zz1;
            const RealType k_p01 = (k_pp1 + k_pm1) + k_pz1;
            const RealType k_m02 = (k_mp2 + k_mm2) + k_mz2;
            const RealType k_z02 = (k_zp2 + k_zm2) + k_zz2;
            const RealType k_p02 = (k_pp2 + k_pm2) + k_pz2;

            const RealType k_m10 = (k_mp0 - k_mm0) - uy * k_m00;
            const RealType k_z10 = (k_zp0 - k_zm0) - uy * k_z00;
            const RealType k_p10 = (k_pp0 - k_pm0) - uy * k_p00;
            const RealType k_m11 = (k_mp1 - k_mm1) - uy * k_m01;
            const RealType k_z11 = (k_zp1 - k_zm1) - uy * k_z01;
            const RealType k_p11 = (k_pp1 - k_pm1) - uy * k_p01;
            const RealType k_m12 = (k_mp2 - k_mm2) - uy * k_m02;
            const RealType k_z12 = (k_zp2 - k_zm2) - uy * k_z02;
            const RealType k_p12 = (k_pp2 - k_pm2) - uy * k_p02;

            const RealType k_m20 = (k_mp0 + k_mm0) - 2.f * uy * (k_mp0 - k_mm0) + uy * uy * k_m00;
            const RealType k_z20 = (k_zp0 + k_zm0) - 2.f * uy * (k_zp0 - k_zm0) + uy * uy * k_z00;
            const RealType k_p20 = (k_pp0 + k_pm0) - 2.f * uy * (k_pp0 - k_pm0) + uy * uy * k_p00;
            const RealType k_m21 = (k_mp1 + k_mm1) - 2.f * uy * (k_mp1 - k_mm1) + uy * uy * k_m01;
            const RealType k_z21 = (k_zp1 + k_zm1) - 2.f * uy * (k_zp1 - k_zm1) + uy * uy * k_z01;
            const RealType k_p21 = (k_pp1 + k_pm1) - 2.f * uy * (k_pp1 - k_pm1) + uy * uy * k_p01;
            const RealType k_m22 = (k_mp2 + k_mm2) - 2.f * uy * (k_mp2 - k_mm2) + uy * uy * k_m02;
            const RealType k_z22 = (k_zp2 + k_zm2) - 2.f * uy * (k_zp2 - k_zm2) + uy * uy * k_z02;
            const RealType k_p22 = (k_pp2 + k_pm2) - 2.f * uy * (k_pp2 - k_pm2) + uy * uy * k_p02;

            //Eq G2015(45)
            // third part of the central moments transformation
            const RealType k_000 = (k_p00 + k_m00) + k_z00;
            const RealType k_001 = (k_p01 + k_m01) + k_z01;
            const RealType k_002 = (k_p02 + k_m02) + k_z02;
            const RealType k_010 = (k_p10 + k_m10) + k_z10;
            const RealType k_011 = (k_p11 + k_m11) + k_z11;
            const RealType k_012 = (k_p12 + k_m12) + k_z12;
            const RealType k_020 = (k_p20 + k_m20) + k_z20;
            const RealType k_021 = (k_p21 + k_m21) + k_z21;
            const RealType k_022 = (k_p22 + k_m22) + k_z22;

            const RealType k_100 = (k_p00 - k_m00) - ux * k_000;
            const RealType k_101 = (k_p01 - k_m01) - ux * k_001;
            const RealType k_102 = (k_p02 - k_m02) - ux * k_002;
            const RealType k_110 = (k_p10 - k_m10) - ux * k_010;
            const RealType k_111 = (k_p11 - k_m11) - ux * k_011;
            const RealType k_112 = (k_p12 - k_m12) - ux * k_012;
            const RealType k_120 = (k_p20 - k_m20) - ux * k_020;
            const RealType k_121 = (k_p21 - k_m21) - ux * k_021;
            const RealType k_122 = (k_p22 - k_m22) - ux * k_022;

            const RealType k_200 = (k_p00 + k_m00) - 2.f * ux * (k_p00 - k_m00) + ux * ux * k_000;
            const RealType k_201 = (k_p01 + k_m01) - 2.f * ux * (k_p01 - k_m01) + ux * ux * k_001;
            const RealType k_202 = (k_p02 + k_m02) - 2.f * ux * (k_p02 - k_m02) + ux * ux * k_002;
            const RealType k_210 = (k_p10 + k_m10) - 2.f * ux * (k_p10 - k_m10) + ux * ux * k_010;
            const RealType k_211 = (k_p11 + k_m11) - 2.f * ux * (k_p11 - k_m11) + ux * ux * k_011;
            const RealType k_212 = (k_p12 + k_m12) - 2.f * ux * (k_p12 - k_m12) + ux * ux * k_012;
            const RealType k_220 = (k_p20 + k_m20) - 2.f * ux * (k_p20 - k_m20) + ux * ux * k_020;
            const RealType k_221 = (k_p21 + k_m21) - 2.f * ux * (k_p21 - k_m21) + ux * ux * k_021;
            const RealType k_222 = (k_p22 + k_m22) - 2.f * ux * (k_p22 - k_m22) + ux * ux * k_022;

            //------------------------------------------------------------------------------------
            //------------------------------ CENTRAL MOM. TO CUMULANTS ---------------------------
            //------------------------------------------------------------------------------------


            //-------------Until order 3 proportional to the central moments---------------

            //Eq G2015(47)
            const RealType C_110 = k_110;
            const RealType C_101 = k_101;
            const RealType C_011 = k_011;

            //Eq G2015(48)
            const RealType C_200 = k_200;
            const RealType C_020 = k_020;
            const RealType C_002 = k_002;

            //Eq G2015(49)
            const RealType C_120 = k_120;
            const RealType C_012 = k_012;
            const RealType C_201 = k_201;

            const RealType C_210 = k_210;
            const RealType C_021 = k_021;
            const RealType C_102 = k_102;

            //Eq G2015(50)
            const RealType C_111 = k_111;

            //Eq G2015(51)
            const RealType C_211 = k_211 - (k_200 * k_011 + 2.f * k_101 * k_110) / rho;
            const RealType C_121 = k_121 - (k_020 * k_101 + 2.f * k_110 * k_011) / rho;
            const RealType C_112 = k_112 - (k_002 * k_110 + 2.f * k_011 * k_101) / rho;

            //Eq G2015(52)
            const RealType C_220 = k_220 - (k_020 * k_200 + 2.f * k_110 * k_110) / rho;
            const RealType C_022 = k_022 - (k_002 * k_020 + 2.f * k_011 * k_011) / rho;
            const RealType C_202 = k_202 - (k_200 * k_002 + 2.f * k_101 * k_101) / rho;

            //Eq G2015(53)
            const RealType C_122 = k_122 - (k_020 * k_102 + k_002 * k_120 + 4.f * k_011 * k_111 +
                                            2.f * (k_110 * k_012 + k_101 * k_021)) / rho;
            const RealType C_212 = k_212 - (k_002 * k_210 + k_200 * k_012 + 4.f * k_101 * k_111 +
                                            2.f * (k_011 * k_201 + k_110 * k_102)) / rho;
            const RealType C_221 = k_221 - (k_200 * k_021 + k_020 * k_201 + 4.f * k_110 * k_111 +
                                            2.f * (k_101 * k_120 + k_011 * k_210)) / rho;
            //Eq G2015(54)
            const RealType C_222 = k_222 - (4.f * k_111 * k_111 + k_200 * k_022 + k_020 * k_202 + k_002 * k_220 +
                                            4.f * (k_011 * k_211 + k_101 * k_121 + k_110 * k_112) +
                                            2.f * (k_120 * k_102 + k_210 * k_012 + k_201 * k_021)) / rho +
                                   (16.f * k_110 * k_101 * k_011 +
                                    4.f * (k_101 * k_101 * k_020 + k_011 * k_011 * k_200 + k_110 * k_110 * k_002) +
                                    2.f * k_200 * k_020 * k_002) / rho / rho;

            //------------------------------------------------------------------------------------
            // -------------------------------------COLLISION-------------------------------------
            //------------------------------------------------------------------------------------

            //  RELAX RATE G2015(103) //2017 diff

            const RealType omega1 = omg; //1.f / (3.f * 1.f / omg + 0.5f);
            const RealType omega2 = 1.f;
            const RealType omega3 = 1.f;
            const RealType omega4 = 1.f;
            const RealType omega5 = 1.f;
            const RealType omega6 = 1.f;
            const RealType omega7 = 1.f;
            const RealType omega8 = 1.f;
            const RealType omega9 = 1.f;
            const RealType omega10 = 1.f;

            //Eq G2015(58)
            const RealType Dxu = -omega1 * 0.5f / rho * (2.f * C_200 - C_020 - C_002) -
                                 omega2 * 0.5f / rho * (C_200 + C_020 + C_002 - k_000); // -(-1-rho));

            //Eq G2015(59)
            const RealType Dyv = Dxu + 3.f * omega1 * 0.5f / rho * (C_200 - C_020);

            //Eq G2015(60)
            const RealType Dzw = Dxu + 3.f * omega1 * 0.5f / rho * (C_200 - C_002);

            //------------------------------------------------------------------------------------

            //Eq G2015(55)
            const RealType Cs_110 = (1.f - omega1) * C_110;
            //Eq G2015(56)
            const RealType Cs_101 = (1.f - omega1) * C_101;
            //Eq G2015(57)
            const RealType Cs_011 = (1.f - omega1) * C_011;

            //---------------------------------------------------------------------------------

            //Eq G2015(61, 62, 63)
            const RealType Eq61RHS = (1.f - omega1) * (C_200 - C_020) -
                                     3.f * rho * (1.f - omega1 * 0.5f) * (ux * ux * Dxu - uy * uy * Dyv);
            const RealType Eq64RHS = (1.f - omega1) * (C_200 - C_002) -
                                     3.f * rho * (1.f - omega1 * 0.5f) * (ux * ux * Dxu - uz * uz * Dzw);
            const RealType Eq65RHS = k_000 * omega2 + (1.f - omega2) * (C_200 + C_020 + C_002) -
                                     3.f * rho * (1.f - omega2 / 2.f) * (ux * ux * Dxu + uy * uy * Dyv + uz * uz * Dzw);

            const RealType Cs_200 = 1.f/3.f * (Eq61RHS + Eq64RHS + Eq65RHS);
            const RealType Cs_020 = 1.f/3.f * (Eq64RHS - 2.f * Eq61RHS + Eq65RHS);
            const RealType Cs_002 = 1.f/3.f * (Eq61RHS - 2.f * Eq64RHS + Eq65RHS);

            //Eq G2015(64, 67)
            const RealType Cs_120 = (-C_102 - C_120) * omega3 * 0.5f + (C_102 - C_120) * omega4 * 0.5f + C_120;
            const RealType Cs_102 = (-C_102 - C_120) * omega3 * 0.5f + (-C_102 + C_120) * omega4 * 0.5f + C_102;

            //Eq G2015(65, 68)
            const RealType Cs_210 = (-C_012 - C_210) * omega3 * 0.5f + (C_012 - C_210) * omega4 * 0.5f + C_210;
            const RealType Cs_012 = (-C_012 - C_210) * omega3 * 0.5f + (-C_012 + C_210) * omega4 * 0.5f + C_012;

            //Eq G2015(66, 69)
            const RealType Cs_021 = (-C_021 - C_201) * omega3 * 0.5f + (-C_021 + C_201) * omega4 * 0.5f + C_021;
            const RealType Cs_201 = (-C_021 - C_201) * omega3 * 0.5f + (C_021 - C_201) * omega4 * 0.5f + C_201;

            //Eq G2015(70)
            const RealType Cs_111 = (1.f - omega5) * C_111;

            //Eq G2015(71, 72, 73) //2017 diff
            const RealType Eq71RHS =  (1.f-omega6)*(C_220-2.f*C_202+C_022);
            const RealType Eq72RHS =  (1.f-omega6)*(C_220+C_202-2.f*C_022);
            const RealType Eq73RHS =  (1.f-omega7)*(C_220+C_202+C_022);
            
            const RealType Cs_220 = 1.f/3.f*(Eq71RHS+Eq72RHS+Eq73RHS);
            const RealType Cs_202 = 1.f/3.f*(-Eq71RHS+Eq73RHS);
            const RealType Cs_022 = 1.f/3.f*(-Eq72RHS+Eq73RHS);
            
            // Eq 46-48 //2017 diff
            const RealType Cs_211 =  (1.f-omega8)*C_211;
            const RealType Cs_121 =  (1.f-omega8)*C_121;
            const RealType Cs_112 =  (1.f-omega8)*C_112;
            // Eqs 49-52
            const RealType Cs_221 = (1.f-omega9)*C_221;
            const RealType Cs_212 = (1.f-omega9)*C_212;
            const RealType Cs_122 = (1.f-omega9)*C_122;
            const RealType Cs_222 = (1.f-omega10)*C_222;

            //------------------------------------------------------------------------------------
            //------------------------------ CUMULANTS TO CENTRAL MOM. ---------------------------
            //------------------------------------------------------------------------------------

            const RealType ks_000 =  k_000;

            // Permutation again

            //Eq G2015(47) backwards
            const RealType ks_110 = Cs_110 ;
            const RealType ks_101 = Cs_101 ;
            const RealType ks_011 = Cs_011 ;

            //Eq G2015(48) backwards
            const RealType ks_200 = Cs_200 ;
            const RealType ks_020 = Cs_020 ;
            const RealType ks_002 = Cs_002 ;

            //Eq G2015(49) backwards
            const RealType ks_120 = Cs_120 ;
            const RealType ks_012 = Cs_012 ;
            const RealType ks_201 = Cs_201 ;

            const RealType ks_210 = Cs_210;
            const RealType ks_021 = Cs_021;
            const RealType ks_102 = Cs_102;

            //Eq G2015(50) backwards
            const RealType ks_111 = Cs_111 ;

            //Eq. G2015(85, 86, 87)
            const RealType ks_100 =  -k_100;
            const RealType ks_010 =  -k_010;
            const RealType ks_001 =  -k_001;

            //Eq. G2015(81)
            const RealType ks_211 = Cs_211+ (ks_200*ks_011 + 2.f*ks_101*ks_110)/rho;
            const RealType ks_121 = Cs_121+ (ks_020*ks_101 + 2.f*ks_110*ks_011)/rho;
            const RealType ks_112 = Cs_112+ (ks_002*ks_110 + 2.f*ks_011*ks_101)/rho;

            //Eq. G2015(82)
            const RealType ks_220 = Cs_220+ (ks_020*ks_200 + 2.f*ks_110*ks_110)/rho;
            const RealType ks_022 = Cs_022+ (ks_002*ks_020 + 2.f*ks_011*ks_011)/rho;
            const RealType ks_202 = Cs_202+ (ks_200*ks_002 + 2.f*ks_101*ks_101)/rho;

            //Eq. G2015(83)
            const RealType ks_122 = Cs_122+ (ks_020*ks_102 + ks_002*ks_120+ 4.f*ks_011*ks_111+ 2.f*(ks_110*ks_012+ks_101*ks_021))/rho;
            const RealType ks_212 = Cs_212+ (ks_002*ks_210 + ks_200*ks_012+ 4.f*ks_101*ks_111+ 2.f*(ks_011*ks_201+ks_110*ks_102))/rho;
            const RealType ks_221 = Cs_221+ (ks_200*ks_021 + ks_020*ks_201+ 4.f*ks_110*ks_111+ 2.f*(ks_101*ks_120+ks_011*ks_210))/rho;
            // gen2nowell.php END

            // Eq. G2015(84)
            const RealType ks_222 = Cs_222 + (4.f*ks_111*ks_111 + ks_200*ks_022 + ks_020*ks_202 + ks_002*ks_220 + 4.f*(ks_011*ks_211 + ks_101*ks_121 + ks_110*ks_112)
                                           + 2.f*(ks_120*ks_102 + ks_210*ks_012 + ks_201*ks_021))/rho
                                 - (16.f*ks_110*ks_101*ks_011 + 4.f*(ks_101*ks_101*ks_020 + ks_011*ks_011*ks_200 + ks_110*ks_110*ks_002) + 2.f*ks_200*ks_020*ks_002)/rho/rho;



            //------------------------------------------------------------------------------------
            //----------------------- TRANSFORM TO DISTRIBUTION FUNCTION -------------------------
            //------------------------------------------------------------------------------------

            //Eq G2015(88)
            const RealType ks_z00 = ks_000 * (1.f - ux * ux) - 2.f * ux * ks_100 - ks_200;
            const RealType ks_z01 = ks_001 * (1.f - ux * ux) - 2.f * ux * ks_101 - ks_201;
            const RealType ks_z02 = ks_002 * (1.f - ux * ux) - 2.f * ux * ks_102 - ks_202;
            const RealType ks_z10 = ks_010 * (1.f - ux * ux) - 2.f * ux * ks_110 - ks_210;
            const RealType ks_z11 = ks_011 * (1.f - ux * ux) - 2.f * ux * ks_111 - ks_211;
            const RealType ks_z12 = ks_012 * (1.f - ux * ux) - 2.f * ux * ks_112 - ks_212;
            const RealType ks_z20 = ks_020 * (1.f - ux * ux) - 2.f * ux * ks_120 - ks_220;
            const RealType ks_z21 = ks_021 * (1.f - ux * ux) - 2.f * ux * ks_121 - ks_221;
            const RealType ks_z22 = ks_022 * (1.f - ux * ux) - 2.f * ux * ks_122 - ks_222;

            //Eq  G2015(89)
            const RealType ks_m00 = (ks_000 * (ux * ux - ux) + ks_100 * (2.f * ux - 1.f) + ks_200) * 0.5f;
            const RealType ks_m01 = (ks_001 * (ux * ux - ux) + ks_101 * (2.f * ux - 1.f) + ks_201) * 0.5f;
            const RealType ks_m02 = (ks_002 * (ux * ux - ux) + ks_102 * (2.f * ux - 1.f) + ks_202) * 0.5f;
            const RealType ks_m10 = (ks_010 * (ux * ux - ux) + ks_110 * (2.f * ux - 1.f) + ks_210) * 0.5f;
            const RealType ks_m11 = (ks_011 * (ux * ux - ux) + ks_111 * (2.f * ux - 1.f) + ks_211) * 0.5f;
            const RealType ks_m12 = (ks_012 * (ux * ux - ux) + ks_112 * (2.f * ux - 1.f) + ks_212) * 0.5f;
            const RealType ks_m20 = (ks_020 * (ux * ux - ux) + ks_120 * (2.f * ux - 1.f) + ks_220) * 0.5f;
            const RealType ks_m21 = (ks_021 * (ux * ux - ux) + ks_121 * (2.f * ux - 1.f) + ks_221) * 0.5f;
            const RealType ks_m22 = (ks_022 * (ux * ux - ux) + ks_122 * (2.f * ux - 1.f) + ks_222) * 0.5f;

            //Eq  G2015(90)
            const RealType ks_p00 = (ks_000 * (ux * ux + ux) + ks_100 * (2.f * ux + 1.f) + ks_200) * 0.5f;
            const RealType ks_p01 = (ks_001 * (ux * ux + ux) + ks_101 * (2.f * ux + 1.f) + ks_201) * 0.5f;
            const RealType ks_p02 = (ks_002 * (ux * ux + ux) + ks_102 * (2.f * ux + 1.f) + ks_202) * 0.5f;
            const RealType ks_p10 = (ks_010 * (ux * ux + ux) + ks_110 * (2.f * ux + 1.f) + ks_210) * 0.5f;
            const RealType ks_p11 = (ks_011 * (ux * ux + ux) + ks_111 * (2.f * ux + 1.f) + ks_211) * 0.5f;
            const RealType ks_p12 = (ks_012 * (ux * ux + ux) + ks_112 * (2.f * ux + 1.f) + ks_212) * 0.5f;
            const RealType ks_p20 = (ks_020 * (ux * ux + ux) + ks_120 * (2.f * ux + 1.f) + ks_220) * 0.5f;
            const RealType ks_p21 = (ks_021 * (ux * ux + ux) + ks_121 * (2.f * ux + 1.f) + ks_221) * 0.5f;
            const RealType ks_p22 = (ks_022 * (ux * ux + ux) + ks_122 * (2.f * ux + 1.f) + ks_222) * 0.5f;

            //Eq G2015(91)
            const RealType ks_mz0 = ks_m00 * (1.f - uy * uy) - 2.f * uy * ks_m10 - ks_m20;
            const RealType ks_mz1 = ks_m01 * (1.f - uy * uy) - 2.f * uy * ks_m11 - ks_m21;
            const RealType ks_mz2 = ks_m02 * (1.f - uy * uy) - 2.f * uy * ks_m12 - ks_m22;
            const RealType ks_zz0 = ks_z00 * (1.f - uy * uy) - 2.f * uy * ks_z10 - ks_z20;
            const RealType ks_zz1 = ks_z01 * (1.f - uy * uy) - 2.f * uy * ks_z11 - ks_z21;
            const RealType ks_zz2 = ks_z02 * (1.f - uy * uy) - 2.f * uy * ks_z12 - ks_z22;
            const RealType ks_pz0 = ks_p00 * (1.f - uy * uy) - 2.f * uy * ks_p10 - ks_p20;
            const RealType ks_pz1 = ks_p01 * (1.f - uy * uy) - 2.f * uy * ks_p11 - ks_p21;
            const RealType ks_pz2 = ks_p02 * (1.f - uy * uy) - 2.f * uy * ks_p12 - ks_p22;

            //Eq  G2015(92)
            const RealType ks_mm0 = (ks_m00 * (uy * uy - uy) + ks_m10 * (2.f * uy - 1.f) + ks_m20) * 0.5f;
            const RealType ks_mm1 = (ks_m01 * (uy * uy - uy) + ks_m11 * (2.f * uy - 1.f) + ks_m21) * 0.5f;
            const RealType ks_mm2 = (ks_m02 * (uy * uy - uy) + ks_m12 * (2.f * uy - 1.f) + ks_m22) * 0.5f;
            const RealType ks_zm0 = (ks_z00 * (uy * uy - uy) + ks_z10 * (2.f * uy - 1.f) + ks_z20) * 0.5f;
            const RealType ks_zm1 = (ks_z01 * (uy * uy - uy) + ks_z11 * (2.f * uy - 1.f) + ks_z21) * 0.5f;
            const RealType ks_zm2 = (ks_z02 * (uy * uy - uy) + ks_z12 * (2.f * uy - 1.f) + ks_z22) * 0.5f;
            const RealType ks_pm0 = (ks_p00 * (uy * uy - uy) + ks_p10 * (2.f * uy - 1.f) + ks_p20) * 0.5f;
            const RealType ks_pm1 = (ks_p01 * (uy * uy - uy) + ks_p11 * (2.f * uy - 1.f) + ks_p21) * 0.5f;
            const RealType ks_pm2 = (ks_p02 * (uy * uy - uy) + ks_p12 * (2.f * uy - 1.f) + ks_p22) * 0.5f;

            //Eq G2015(93)
            const RealType ks_mp0 = (ks_m00 * (uy * uy + uy) + ks_m10 * (2.f * uy + 1.f) + ks_m20) * 0.5f;
            const RealType ks_mp1 = (ks_m01 * (uy * uy + uy) + ks_m11 * (2.f * uy + 1.f) + ks_m21) * 0.5f;
            const RealType ks_mp2 = (ks_m02 * (uy * uy + uy) + ks_m12 * (2.f * uy + 1.f) + ks_m22) * 0.5f;
            const RealType ks_zp0 = (ks_z00 * (uy * uy + uy) + ks_z10 * (2.f * uy + 1.f) + ks_z20) * 0.5f;
            const RealType ks_zp1 = (ks_z01 * (uy * uy + uy) + ks_z11 * (2.f * uy + 1.f) + ks_z21) * 0.5f;
            const RealType ks_zp2 = (ks_z02 * (uy * uy + uy) + ks_z12 * (2.f * uy + 1.f) + ks_z22) * 0.5f;
            const RealType ks_pp0 = (ks_p00 * (uy * uy + uy) + ks_p10 * (2.f * uy + 1.f) + ks_p20) * 0.5f;
            const RealType ks_pp1 = (ks_p01 * (uy * uy + uy) + ks_p11 * (2.f * uy + 1.f) + ks_p21) * 0.5f;
            const RealType ks_pp2 = (ks_p02 * (uy * uy + uy) + ks_p12 * (2.f * uy + 1.f) + ks_p22) * 0.5f;

            //Eq G2015(94)
            f_mmz = ks_mm0 * (1.f - uz * uz) - 2.f * uz * ks_mm1 - ks_mm2;
            f_mzz = ks_mz0 * (1.f - uz * uz) - 2.f * uz * ks_mz1 - ks_mz2;
            f_mpz = ks_mp0 * (1.f - uz * uz) - 2.f * uz * ks_mp1 - ks_mp2;
            f_zmz = ks_zm0 * (1.f - uz * uz) - 2.f * uz * ks_zm1 - ks_zm2;
            f_zzz = ks_zz0 * (1.f - uz * uz) - 2.f * uz * ks_zz1 - ks_zz2;
            f_zpz = ks_zp0 * (1.f - uz * uz) - 2.f * uz * ks_zp1 - ks_zp2;
            f_pmz = ks_pm0 * (1.f - uz * uz) - 2.f * uz * ks_pm1 - ks_pm2;
            f_pzz = ks_pz0 * (1.f - uz * uz) - 2.f * uz * ks_pz1 - ks_pz2;
            f_ppz = ks_pp0 * (1.f - uz * uz) - 2.f * uz * ks_pp1 - ks_pp2;

            //Eq  G2015(95)
            f_mmm = (ks_mm0 * (uz * uz - uz) + ks_mm1 * (2.f * uz - 1.f) + ks_mm2) * 0.5f;
            f_mzm = (ks_mz0 * (uz * uz - uz) + ks_mz1 * (2.f * uz - 1.f) + ks_mz2) * 0.5f;
            f_mpm = (ks_mp0 * (uz * uz - uz) + ks_mp1 * (2.f * uz - 1.f) + ks_mp2) * 0.5f;
            f_zmm = (ks_zm0 * (uz * uz - uz) + ks_zm1 * (2.f * uz - 1.f) + ks_zm2) * 0.5f;
            f_zzm = (ks_zz0 * (uz * uz - uz) + ks_zz1 * (2.f * uz - 1.f) + ks_zz2) * 0.5f;
            f_zpm = (ks_zp0 * (uz * uz - uz) + ks_zp1 * (2.f * uz - 1.f) + ks_zp2) * 0.5f;
            f_pmm = (ks_pm0 * (uz * uz - uz) + ks_pm1 * (2.f * uz - 1.f) + ks_pm2) * 0.5f;
            f_pzm = (ks_pz0 * (uz * uz - uz) + ks_pz1 * (2.f * uz - 1.f) + ks_pz2) * 0.5f;
            f_ppm = (ks_pp0 * (uz * uz - uz) + ks_pp1 * (2.f * uz - 1.f) + ks_pp2) * 0.5f;

            //Eq  G2015(96)
            f_mmp = (ks_mm0 * (uz * uz + uz) + ks_mm1 * (2.f * uz + 1.f) + ks_mm2) * 0.5f;
            f_mzp = (ks_mz0 * (uz * uz + uz) + ks_mz1 * (2.f * uz + 1.f) + ks_mz2) * 0.5f;
            f_mpp = (ks_mp0 * (uz * uz + uz) + ks_mp1 * (2.f * uz + 1.f) + ks_mp2) * 0.5f;
            f_zmp = (ks_zm0 * (uz * uz + uz) + ks_zm1 * (2.f * uz + 1.f) + ks_zm2) * 0.5f;
            f_zzp = (ks_zz0 * (uz * uz + uz) + ks_zz1 * (2.f * uz + 1.f) + ks_zz2) * 0.5f;
            f_zpp = (ks_zp0 * (uz * uz + uz) + ks_zp1 * (2.f * uz + 1.f) + ks_zp2) * 0.5f;
            f_pmp = (ks_pm0 * (uz * uz + uz) + ks_pm1 * (2.f * uz + 1.f) + ks_pm2) * 0.5f;
            f_pzp = (ks_pz0 * (uz * uz + uz) + ks_pz1 * (2.f * uz + 1.f) + ks_pz2) * 0.5f;
            f_ppp = (ks_pp0 * (uz * uz + uz) + ks_pp1 * (2.f * uz + 1.f) + ks_pp2) * 0.5f;

            // back to df_post
            df_post_view(i, j, k, 0) = f_zzz; // 0
            df_post_view(i, j, k, 1) = f_pzz; // 1
            df_post_view(i, j, k, 2) = f_mzz; // 2
            df_post_view(i, j, k, 6) = f_zpz; // 6
            df_post_view(i, j, k, 5) = f_zmz; // 5
            df_post_view(i, j, k, 4) = f_zzp; // 4
            df_post_view(i, j, k, 3) = f_zzm; // 3

            df_post_view(i, j, k, 12) = f_ppz; // 12
            df_post_view(i, j, k, 11) = f_mmz; // 11
            df_post_view(i, j, k, 16) = f_pmz; // 16
            df_post_view(i, j, k, 15) = f_mpz; // 15
            df_post_view(i, j, k, 9) = f_pzp; // 9
            df_post_view(i, j, k, 10) = f_mzm; // 10
            df_post_view(i, j, k, 7) = f_pzm; // 7
            df_post_view(i, j, k, 8) = f_mzp; // 8
            df_post_view(i, j, k, 17) = f_zpp; // 17
            df_post_view(i, j, k, 18) = f_zmm; // 18
            df_post_view(i, j, k, 13) = f_zpm; // 13
            df_post_view(i, j, k, 14) = f_zmp; // 14

            df_post_view(i, j, k, 26) = f_ppp; // 26
            df_post_view(i, j, k, 25) = f_mmm; // 25
            df_post_view(i, j, k, 22) = f_ppm; // 22
            df_post_view(i, j, k, 21) = f_mmp; // 21
            df_post_view(i, j, k, 20) = f_pmp; // 20
            df_post_view(i, j, k, 19) = f_mpm; // 19
            df_post_view(i, j, k, 23) = f_pmm; // 23
            df_post_view(i, j, k, 24) = f_mpp; // 24

        };


        auto coll = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<3, int> &i )mutable
        {

            if (mesh_view(i.x(), i.y(), i.z()) != 0) {

                cum(i.x(), i.y(), i.z());
            }

        };


        TNL::Containers::StaticArray<3, int> begin{0, 0, 0};
        TNL::Containers::StaticArray<3, int> end{Constants->dimX_int, Constants->dimY_int, Constants->dimZ_int};
        parallelFor<DeviceType>(begin, end, coll);
    }


};


#endif
