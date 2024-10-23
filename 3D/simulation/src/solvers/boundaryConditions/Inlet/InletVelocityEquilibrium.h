//
// Created by stloufra on 10/30/23.
//

#ifndef INLETVELOCITYEQUILIBRIUM_H
#define INLETVELOCITYEQUILIBRIUM_H

#include "../../../traits/LBMTraits.h"

template<typename MODELDATA>
struct InletVelocityEquilibrium {
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;

    static void inlet(LBMDataPointer &Data, LBMConstantsPointer &Constants) {


        auto inlet_view = Data->meshBoundaryInlet.getView();

        auto rho_view = Data->rho.getView();
        auto u_view = Data->u.getView();

        auto df_view = Data->df.getView();
        auto df_post_view = Data->df_post.getView();

        const auto Nvel = Constants->Nvel;
        auto ny = Constants->ny;

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

        auto bb_inlet = [=]
        __cuda_callable__(
        const TNL::Containers::StaticArray<1, int> &nod  ) mutable
        {
            Vertex vert = inlet_view[nod.x()].vertex;
            Vector norm = inlet_view[nod.x()].normal;
            Vector velc = inlet_view[nod.x()].velocity;

            int i=vert.x;
            int j=vert.y;
            int k=vert.z;

            const RealType uxb = velc.x();
            const RealType uyb = velc.y();
            const RealType uzb = velc.z();

            const RealType rho = rho_view(i, j, k);

            Vector xp( 1.0f, 0.0f, 0.0f ), xm( -1.0f, 0.0f, 0.0f );

             if ( norm == xp){
                 RealType ux1 = u_view(i - 1 , j, k).x();
                 RealType ux2 = u_view(i - 2 , j, k).x();

                 RealType duxb = (uxb-ux1); //dx = 1
                 RealType dux1 = (ux1-ux2);

                 RealType rhob = rho_view(i-1,j,k) + ny/3*(dux1 - duxb); //1/cs^2 = 3

                 for(int v=0; v<Nvel; v++){
                     df_view(i,j,k,v)= f_equilibrium_defined(uxb, uyb, uzb, rhob, v);
                 }
             }
             if ( norm == xm){
                 RealType ux1 = u_view(i + 1 , j, k).x();
                 RealType ux2 = u_view(i + 2 , j, k).x();

                 RealType duxb = (ux1-uxb); //dx = 1
                 RealType dux1 = (ux2-ux1);

                 RealType rhob = rho_view(i+1,j,k) + ny/3*(duxb - dux1); //1/cs^2 = 3

                 for(int v=0; v<Nvel; v++){
                     df_view(i,j,k,v)= f_equilibrium_defined(uxb, uyb, uzb, rhob, v);
                 }

             }
             else{
               printf("Not yet supported inlet.\n");
             }


        };


        parallelFor<DeviceType>(0, Constants->inlet_num, bb_inlet);

    }

};

#endif
