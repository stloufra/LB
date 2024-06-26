#ifndef SOLVER_LES_H
#define SOLVER_LES_H

#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>
#include <omp.h>
#include "Mesher.h"
#include <TNL/Containers/NDArray.h>
#include <TNL/Containers/StaticArray.h>
#include <TNL/Containers/StaticVector.h>
#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Algorithms/reduce.h>


using namespace TNL;
using namespace TNL::Containers;
using namespace TNL::Algorithms;

using DeviceTypeHost = TNL::Devices::Host;


template< typename RealType, typename DeviceType >
class Solver_LES
{

using RealType1 = float;

using ArrayType2D =  TNL::Containers::NDArray<  RealType,
                                                SizesHolder<int, 0, 0 >,
                                                std::index_sequence< 0, 1 > ,
                                                DeviceType>;

using ArrayType2D_Host =  TNL::Containers::NDArray<  RealType,
                                                SizesHolder<int, 0, 0 >,
                                                std::index_sequence< 0, 1 > ,
                                                DeviceTypeHost>;                                              

using ArrayType3D =  TNL::Containers::NDArray<  RealType,
                                                SizesHolder<int, 0, 0, 0 >,
                                                std::index_sequence< 0, 1, 2> ,
                                                DeviceType>;

using TensorType = TNL::Containers::StaticVector<2, TNL::Containers::StaticVector< 2, RealType1 >>;


public:
  Solver_LES()=delete;

  Solver_LES(int Ny_in, int Nx_in, Mesher<RealType, DeviceTypeHost> mesh_in):
  rho(),
  ux(),
  uy(),
  omega(),
  ux0er(),
  uy0er(),
  df(),
  df_post(),
  mesh(),
  velocities_x(),
  velocities_y()
  {   
    Nx = Nx_in;
    Ny = Ny_in;

    rho.setSizes(Ny, Nx);
    ux.setSizes(Ny, Nx);
    uy.setSizes(Ny, Nx);
    omega.setSizes(Ny, Nx);
    ux0er.setSizes(Ny, Nx);
    uy0er.setSizes(Ny, Nx);
    df.setSizes(Ny, Nx, Nvel);
    df_post.setSizes(Ny, Nx, Nvel);
              
    mesh = mesh_in.mesh;
    velocities_x = mesh_in.velocities_x;
    velocities_y = mesh_in.velocities_y;
  }

  ~Solver_LES() {}

  void initialization_eq(RealType rho0, RealType ux0, RealType uy0, RealType Fx0, RealType Fy0, bool gravity)
  {

    auto rho_view = rho.getView();
    auto ux_view = ux.getView();
    auto uy_view = uy.getView();
    auto omega_view = omega.getView();
    auto mesh_view = mesh.getView();
    auto df_view = df.getView();
    auto df_post_view = df_post.getView();


    const RealType Cl = this->Cl;
    const RealType Ct = this->Ct;
    const RealType Cm = this->Cm;
    const RealType Cu_inverse = this->Cu_inverse;
    const RealType Cu = this->Cu;
    const int Nvel = this-> Nvel;
    const int (&cx_pos)[9] = this->cx_pos;
    const int (&cy_pos)[9] = this->cy_pos;
    const RealType (&weight)[9] = this->weight;
    RealType omega_0 = this->omega_0;
    

    auto f_equilibrium = [=] __cuda_callable__ ( const int &j,const int &i, const int &k ) mutable 
    {
      RealType cu, u2;

      cu = cx_pos[k]*ux_view(j,i) + cy_pos[k]*uy_view(j,i);
      u2 = ux_view(j,i)*ux_view(j,i) + uy_view(j,i)*uy_view(j,i);

      return weight[k]*rho_view(j,i)*(1.f + 3.f*cu + 4.5f*cu*cu - 1.5f*u2);
    };

    auto init_variables = [=] __cuda_callable__ ( const StaticArray< 2, int >& i  ) mutable 
    {
      if( mesh_view(i.y()+1, i.x()+1) !=0 )
      {
        rho_view(i.y(),i.x())=rho0/Cm*Cl*Cl*Cl;
        ux_view(i.y(),i.x()) = ux0*Cu_inverse;
        uy_view(i.y(),i.x()) = uy0*Cu_inverse;
        omega_view(i.y(),i.x()) = omega_0;
      }

      else if(mesh_view(i.y()+1, i.x()+1)==0)
      {
        rho_view(i.y(),i.x())= 0.f;
        ux_view(i.y(),i.x()) = 0.f;
        uy_view(i.y(),i.x()) = 0.f;
      }
    };

    auto init_df = [=] __cuda_callable__ ( const StaticArray< 3, int >& i  ) mutable 
    {
      if( mesh_view(i.y()+1,i.x()+1) !=0 )
      {
        df_view(i.y(),i.x(), i.z()) = f_equilibrium(i.y(),i.x(), i.z());
        df_post_view(i.y(),i.x(), i.z()) = f_equilibrium(i.y(),i.x(), i.z());

      }
      else if(mesh_view(i.y()+1,i.x()+1)==0)
      {
        df_view(i.y(),i.x(), i.z()) = 0.f;
        df_post_view(i.y(),i.x(), i.z()) = 0.f;
      }
    };

    
    StaticArray< 2, int > begin1{ 0, 0};
    StaticArray< 2, int > end1{ Nx, Ny};
    //std::cout << "size: " << end1 << " array size: " << rho_view.getSizes() << " " << mesh_view.getSizes() << std::endl; 
    parallelFor< DeviceType >( begin1, end1, init_variables );
    //printf( "init_variabes \n");

    StaticArray< 3, int > begin2{ 0, 0, 0};
    StaticArray< 3, int > end2{ Nx, Ny, Nvel};
    parallelFor< DeviceType >( begin2, end2, init_df );
    //printf( "init_df \n");


    Fx = Fx0/Cm*Cl*Cl*Ct*Ct;
    Fy = Fy0/Cm*Cl*Cl*Ct*Ct;

    if (gravity == 1)
    {
      g=9.81f/Cl*Ct*Ct;
    }
    else 
    {
      g=0.f;
    }


    std::cout<<"Initialization undergone successfully.\n";

  }

  void collision()
  {
    auto rho_view = rho.getView();
    auto ux_view = ux.getView();
    auto uy_view = uy.getView();
    auto mesh_view = mesh.getView();
    auto df_view = df.getView();
    auto omega_view = omega.getView();
    auto df_post_view = df_post.getView();

    const int Nvel = this-> Nvel;
    const int (&cx_pos)[9] = this->cx_pos;
    const int (&cy_pos)[9] = this->cy_pos;
    const RealType (&weight)[9] = this->weight;



    auto f_equilibrium = [=] __cuda_callable__ ( const int &j, const int &i, const int &k ) mutable 
    {
      RealType cu, u2;

      cu = cx_pos[k]*ux_view(j,i) + cy_pos[k]*uy_view(j,i);
      u2 = ux_view(j,i)*ux_view(j,i) + uy_view(j,i)*uy_view(j,i);

      return weight[k]*rho_view(j,i)*(1.f + 3.f*cu + 4.5f*cu*cu - 1.5f*u2);
    };

    auto coll = [=] __cuda_callable__ ( const StaticArray< 3, int >& i  ) mutable 
    {
      if (mesh_view(i.y()+1, i.x()+1) == 1)
      {
        RealType feq = f_equilibrium(i.y(), i.x(), i.z());
        df_post_view(i.y(), i.x(), i.z()) = df_view(i.y(),i.x(), i.z()) - (df_view(i.y(), i.x(), i.z())-feq)*omega_view(i.y(), i.x());
      }
    };

    StaticArray< 3, int > begin{ 0, 0, 0};
    StaticArray< 3, int > end{ Nx, Ny, Nvel};
    parallelFor< DeviceType >( begin, end, coll );
    
    //std::cout<<"collision sucesfull.\n";
  }

  void streaming()
  {
    auto df_view = df.getView();
    auto df_post_view = df_post.getView();
    auto mesh_view = mesh.getView();

    const int Nvel = this-> Nvel;
    const int Nx = this-> Nx;
    const int Ny = this-> Ny;
    const int (&cx_pos)[9] = this->cx_pos;
    const int (&cy_pos)[9] = this->cy_pos;
    
    
    auto stream = [=] __cuda_callable__ ( const StaticArray< 3, int >& i  ) mutable 
    {
      if (mesh_view(i.y()+1, i.x()+1) == 1)
      {
        int jd; int id;
        
        jd = i.y() - cy_pos[i.z()];
        id = i.x() - cx_pos[i.z()];

        if(mesh_view(jd+1,id+1)==1 && jd>=0 && jd < Ny && id>=0 && id < Nx)
        {       
          df_view(i.y(), i.x(), i.z())=df_post_view(jd,id,i.z());
        } 
      }

      if (mesh_view(i.y()+1, i.x()+1) == 0)
      {
        df_view(i.y(),i.x(), i.z())=0;
        
      }
    };

    StaticArray< 3, int > begin{ 0, 0, 0};
    StaticArray< 3, int > end{ Nx, Ny, Nvel};
    parallelFor< DeviceType >( begin, end, stream );
  
    //std::cout<<"streaming sucesfull.\n";
    }

  void bounce_back()
  {
    auto df_view = df.getView();
    auto df_post_view = df_post.getView();
    auto mesh_view = mesh.getView();
    auto velocities_x_view = velocities_x.getView();
    auto velocities_y_view = velocities_y.getView();
    auto rho_view = rho.getView();
    auto ux_view = ux.getView();
    auto uy_view = uy.getView();

    const int Nvel = this-> Nvel;
    const int Nx = this-> Nx;
    const int Ny = this-> Ny;
    const int (&cx_pos)[9] = this->cx_pos;
    const int (&cy_pos)[9] = this->cy_pos;
    const RealType (&weight)[9] = this->weight;
    const int (&c_rever)[9] = this->c_rever;

    const RealType Cu_inverse = this->Cu_inverse;


    // MESH - structured bolean values of BC
    // 0 = solid | 1 = fluid | 2 = primitive inlet vertical | 3 = outlet (equilibrium) | 4 = moving wall up | 5 = moving wall down

    //BOUNCE BACK WALL
    auto f_equilibrium = [=] __cuda_callable__ ( const int &j,const int &i, const int &k ) mutable 
    {
      RealType cu, u2;

      cu = cx_pos[k]*ux_view(j,i) + cy_pos[k]*uy_view(j,i);
      u2 = ux_view(j,i)*ux_view(j,i) + uy_view(j,i)*uy_view(j,i);

      return weight[k]*rho_view(j,i)*(1.f + 3.f*cu + 4.5f*cu*cu - 1.5f*u2);
    };

    auto bb_wall = [=] __cuda_callable__ ( const StaticArray< 3, int >& i  ) mutable 
    {
      int dx = cx_pos[i.z()];
      int dy = cy_pos[i.z()];

      if(mesh_view(i.y()+1,i.x()+1) ==1)
      {
        int neighbour_x; int neighbour_y;

        neighbour_x = i.x() - dx;
        neighbour_y = i.y() - dy;

        if(mesh_view(neighbour_y+1,neighbour_x+1)==0)
        {
          df_view(i.y(), i.x(), i.z()) = df_post_view(i.y(), i.x(), c_rever[i.z()]);
        }                 
      }
    };

    // BOUNDARY CONDITIONS

    auto bb_bc = [=] __cuda_callable__ ( const StaticArray< 2, int >& i  ) mutable 
    {
      if(mesh_view(i.y()+1,i.x()+1)==2) //INLET
      {
          if(velocities_x_view(i.y()+1,i.x()+1)>0)
          {
            int x = i.x() + 1;
            df_view(i.y(),x,1)=df_post_view(i.y(),x,3) + 6*weight[1]*rho_view(i.y(),x)*velocities_x_view(i.y()+1,i.x()+1)*Cu_inverse;
            df_view(i.y(),x,5)=df_post_view(i.y(),x,7) + 6*weight[5]*rho_view(i.y(),x)*velocities_x_view(i.y()+1,i.x()+1)*Cu_inverse;
            df_view(i.y(),x,8)=df_post_view(i.y(),x,6) + 6*weight[8]*rho_view(i.y(),x)*velocities_x_view(i.y()+1,i.x()+1)*Cu_inverse;  
          }

          else if(velocities_x_view(i.y()+1,i.x()+1)<0)
          {
            int x = i.x() - 1;
            df_view(i.y(),x,3)=df_post_view(i.y(),x,1) + 6*weight[1]*rho_view(i.y(),x)*velocities_x_view(i.y()+1,i.x()+1)*Cu_inverse;
            df_view(i.y(),x,7)=df_post_view(i.y(),x,5) + 6*weight[5]*rho_view(i.y(),x)*velocities_x_view(i.y()+1,i.x()+1)*Cu_inverse;
            df_view(i.y(),x,6)=df_post_view(i.y(),x,8) + 6*weight[8]*rho_view(i.y(),x)*velocities_x_view(i.y()+1,i.x()+1)*Cu_inverse;  
          }
      }

      else if(mesh_view(i.y()+1,i.x()+1)==3) //OUTLET RIGHT
      {   
        int x = i.x() - 1;
        RealType rho_out = 1.f;

        df_view(i.y(),x,3)=f_equilibrium(i.y(),i.x(), 3);
        df_view(i.y(),x,7)=f_equilibrium(i.y(),i.x(), 7);
        df_view(i.y(),x,6)=f_equilibrium(i.y(),i.x(), 6);  
      } 

      else if(mesh_view(i.y()+1,i.x()+1)==6) //OUTLET LEFT
      {
        int x = i.x() + 1;
        RealType rho_out = 1.f;

        df_view(i.y(),x,5)=f_equilibrium(i.y(),i.x(), 5);
        df_view(i.y(),x,1)=f_equilibrium(i.y(),i.x(), 1);
        df_view(i.y(),x,8)=f_equilibrium(i.y(),i.x(), 8);  
      } 
      //df(i.y(),x,3)=f_equilibrium(rho(i.y(),x), ux(i.y(),x), uy(i.y(),x), 3);
      //df(i.y(),x,7)=f_equilibrium(rho(i.y(),x), ux(i.y(),x), uy(i.y(),x), 7);
      //df(i.y(),x,6)=f_equilibrium(rho(i.y(),x), ux(i.y(),x), uy(i.y(),x), 6);  

      else if(mesh_view(i.y()+1,i.x()+1)==4) //WALL MOVING UP
      {
        if(i.x()>-1 && i.x()<Nx)
        {
          int y = i.y() - 1;
          df_view(y,i.x(),4)=df_post_view(y,i.x(),2);
          df_view(y,i.x(),7)=df_post_view(y,i.x(),5)+6*rho_view(y,i.x())*weight[7]*cx_pos[7]*velocities_x_view(i.y()+1,i.x()+1)*Cu_inverse;
          df_view(y,i.x(),8)=df_post_view(y,i.x(),6)+6*rho_view(y,i.x())*weight[8]*cx_pos[8]*velocities_x_view(i.y()+1,i.x()+1)*Cu_inverse;
        }
      }

      else if(mesh_view(i.y()+1,i.x()+1)==5) //WALL MOVING DOWN
      {   
        if(i.x()>-1 && i.x()<Nx)
        {
          int y = i.y() + 1;
          df_view(y,i.x(),2)=df_post_view(y,i.x(),4);
          df_view(y,i.x(),5)=df_post_view(y,i.x(),7)+6*rho_view(y,i.x())*weight[5]*cx_pos[5]*velocities_x_view(i.y()+1,i.x()+1)*Cu_inverse;
          df_view(y,i.x(),6)=df_post_view(y,i.x(),8)+6*rho_view(y,i.x())*weight[6]*cx_pos[6]*velocities_x_view(i.y()+1,i.x()+1)*Cu_inverse; 
        }
      }
    };


    StaticArray< 3, int > begin1{ 0, 0, 0};
    StaticArray< 3, int > end1{ Nx, Ny, Nvel};
    parallelFor< DeviceType >( begin1, end1, bb_wall);

    StaticArray< 2, int > begin2{ 0, 0};
    StaticArray< 2, int > end2{ Nx, Ny};
    parallelFor< DeviceType >( begin2, end2, bb_bc );
    //std::cout<<"bb sucesfull.\n";
  }

  void postpro()
  {
    auto rho_view = rho.getView();
    auto ux_view = ux.getView();
    auto uy_view = uy.getView();
    auto df_view = df.getView();
    auto mesh_view = mesh.getView();

    RealType tau = this->tau;
    RealType g = this->g;
    RealType Fx = this->Fx;
    RealType Fy = this->Fy;
   
    auto post = [=] __cuda_callable__ ( const StaticArray< 2, int >& i  ) mutable 
    {
      if(mesh_view(i.y()+1,i.x()+1)==1)
      {
        rho_view(i.y(),i.x())=df_view(i.y(),i.x(),0)+df_view(i.y(),i.x(),1)+df_view(i.y(),i.x(),2)+df_view(i.y(),i.x(),3)+df_view(i.y(),i.x(),4)+df_view(i.y(),i.x(),5)+df_view(i.y(),i.x(),6)+df_view(i.y(),i.x(),7)+df_view(i.y(),i.x(),8);
        ux_view(i.y(),i.x())=(df_view(i.y(),i.x(),1)+df_view(i.y(),i.x(),5)+df_view(i.y(),i.x(),8)-df_view(i.y(),i.x(),3)-df_view(i.y(),i.x(),6)-df_view(i.y(),i.x(),7))/rho_view(i.y(),i.x()) + Fx/rho_view(i.y(),i.x())*tau;
        uy_view(i.y(),i.x())=(df_view(i.y(),i.x(),5)+df_view(i.y(),i.x(),6)+df_view(i.y(),i.x(),2)-df_view(i.y(),i.x(),7)-df_view(i.y(),i.x(),8)-df_view(i.y(),i.x(),4))/rho_view(i.y(),i.x()) + Fy/rho_view(i.y(),i.x())*tau - g*tau;      
      }
    };

    StaticArray< 2, int > begin{ 0, 0};
    StaticArray< 2, int > end{ Nx, Ny};
    parallelFor< DeviceType >( begin, end, post );
  }

  void omegaLES(){

      auto rho_view = rho.getView();
      auto ux_view = ux.getView();
      auto uy_view = uy.getView();
      auto omega_view = omega.getView();
      auto mesh_view = mesh.getView();
      auto df_view = df.getView();
      auto df_post_view = df_post.getView();

      const int Nvel = this-> Nvel;
      const int (&cx_pos)[9] = this->cx_pos;
      const int (&cy_pos)[9] = this->cy_pos;
      const int (&c)[9][2] = this->c;
      const RealType (&weight)[9] = this->weight;
      const auto cs4 = this->cs4;

      const auto tau = this->tau;
      const auto Cl = this-> Cl;
      const auto Ct =  this-> Ct;
      const auto CLES = this-> CLES;



      auto f_equilibrium = [=] __cuda_callable__ ( const int &j, const int &i, const int &k ) mutable
      {
          RealType cu, u2;

          cu = cx_pos[k]*ux_view(j,i) + cy_pos[k]*uy_view(j,i);
          u2 = ux_view(j,i)*ux_view(j,i) + uy_view(j,i)*uy_view(j,i);

          return weight[k]*rho_view(j,i)*(1.f + 3.f*cu + 4.5f*cu*cu - 1.5f*u2);
      };



      auto omega = [=]
      __cuda_callable__(
      const TNL::Containers::StaticArray<3, int> &i ) mutable
      {
          if (mesh_view(i.y()+1, i.x()+1) == 1) {
              TensorType pi;

              RealType PI=0.f;

              RealType feq = f_equilibrium(i.y(), i.x(), i.z());

              for (int alpha = 0; alpha < 2; alpha++) {
                  for (int beta = 0; beta < 2; beta++) {
                          pi(alpha)(beta) += c[i.z()][alpha] * c[i.z()][beta] * (df_view(i.y(), i.x(), i.z()) -
                                                                                   feq);

                      PI += pi(alpha)(beta)*pi(alpha)(beta);
                  }
              }

              PI = sqrt(PI);

              RealType tauLES;

              tauLES = tau * 0.5f +
                       sqrt(tau * tau + 2.f * CLES * Cl * Cl / (rho_view(i.y(), i.x()) * cs4 * Ct * Ct) * PI) *
                       0.5f;

              omega_view(i.y(), i.x()) = 1.f / tauLES;
          }
      };

      StaticArray< 3, int > begin{ 0, 0, 0};
      StaticArray< 3, int > end{ Nx, Ny, Nvel};
      parallelFor< DeviceType >( begin, end, omega);

  }

  void Err()
  {
    auto ux_view = ux.getView();
    auto uy_view = uy.getView();
    auto ux0er_view = ux0er.getView();
    auto uy0er_view = uy0er.getView();

    auto ux_ArrayData_view = ux.getStorageArray().getConstView();
    auto uy_ArrayData_view = uy.getStorageArray().getConstView();
    auto ux0er_ArrayData_view = ux0er.getStorageArray().getConstView();
    auto uy0er_ArrayData_view = uy0er.getStorageArray().getConstView();
    
    RealType err1, err2;
    err1=0.0f; err2=0.0f;

    auto fetch1 = [=] __cuda_callable__ ( int i ) ->RealType 
    {
      const RealType& add1 = (ux_ArrayData_view[i]-ux0er_ArrayData_view[i])*(ux_ArrayData_view[i]-ux0er_ArrayData_view[i])+(uy_ArrayData_view[i]-uy0er_ArrayData_view[i])*(uy_ArrayData_view[i]-uy0er_ArrayData_view[i]);
      return sqrt(add1); 
    };

    auto reduction1 = [] __cuda_callable__ ( const RealType& a, const RealType& b )
    {
        return a + b;
    };
    
    err1 = sqrt( reduce<DeviceType, int>( 0,
                                              ux.getStorageArray().getSize(),
                                              fetch1,
                                              reduction1,
                                              0.0 ) );

    auto fetch2 = [=] __cuda_callable__ ( int i ) ->RealType 
    {
      const RealType& add2 = ux_ArrayData_view[i]*ux_ArrayData_view[i] + uy_ArrayData_view[i]*uy_ArrayData_view[i];
      return sqrt(add2); 
    };

    auto reduction2 = [] __cuda_callable__ ( const RealType& a, const RealType& b )
    {
        return a + b;
    };
  
    err2 = sqrt( reduce<DeviceType, int>( 0,
                                              ux.getStorageArray().getSize(),
                                              fetch2,
                                              reduction2,
                                              0.0 ) );
    
    err=err1/err2;
    

    auto copy = [=] __cuda_callable__ ( const StaticArray< 2, int >& i  ) mutable 
    {
        ux0er_view(i.y(),i.x()) = ux_view(i.y(),i.x());
        uy0er_view(i.y(),i.x()) = uy_view(i.y(),i.x());
    };

    StaticArray< 2, int > begin3{ 0, 0};
    StaticArray< 2, int > end3{ Nx, Ny};
    parallelFor< DeviceType >( begin3, end3, copy );
  }

  void output_VTK_lattice()
  {
       
    if(std::is_same_v< DeviceType, TNL::Devices::Cuda > )
    {
      //std::cout<<"\nCuda -> Host output Lattice units\n"<<std::endl;
      ArrayType2D_Host rho_out;
      ArrayType2D_Host ux_out;
      ArrayType2D_Host uy_out;
      
      rho_out = rho;
      ux_out = ux;
      uy_out = uy;

      std::ofstream out_file("results/Lattice_units.vtk");

      out_file << "# vtk DataFile Version 2.0\n" ;
      out_file << "LBE mesh\n" ;
      out_file << "ASCII\n";
      out_file << "DATASET STRUCTURED_POINTS\n" ;
      out_file << "DIMENSIONS "<<Nx<<" "<<Ny<<" 1\n" ;
      out_file << "ASPECT_RATIO 1 1 1\n";
      out_file << "ORIGIN 0 0 0\n";
      out_file << "POINT_DATA "<<Nx*Ny<<"\n";

      out_file << "SCALARS "<<"rho "<<"double 1\n";
      out_file << "LOOKUP_TABLE default\n";

      
      for(int j=0; j < Ny; j++)
      {
          for(int i=0; i < Nx; i++)
          {
              out_file <<rho_out(j,i)<<"\n";
          }
      }
        
      out_file << "VECTORS "<<"U "<<"double\n";

      
      for(int j=0; j < Ny; j++)
      {
          for(int i=0; i <Nx; i++)
          {
              out_file <<ux_out(j,i)<<" "<<uy_out(j,i)<<" 0\n";
          }
      }

      out_file.close();

    } 
    else if (std::is_same<DeviceType, TNL::Devices::Host >::value )
    {
      //std::cout<<"\nHost output Lattice units\n"<<std::endl;
      
      std::ofstream out_file("results/Lattice_units.vtk");

      out_file << "# vtk DataFile Version 2.0\n" ;
      out_file << "LBE mesh\n" ;
      out_file << "ASCII\n";
      out_file << "DATASET STRUCTURED_POINTS\n" ;
      out_file << "DIMENSIONS "<<Nx<<" "<<Ny<<" 1\n" ;
      out_file << "ASPECT_RATIO 1 1 1\n";
      out_file << "ORIGIN 0 0 0\n";
      out_file << "POINT_DATA "<<Nx*Ny<<"\n";

      out_file << "SCALARS "<<"rho "<<"double 1\n";
      out_file << "LOOKUP_TABLE default\n";

      
      for(int j=0; j < Ny; j++)
      {
          for(int i=0; i < Nx; i++)
          {
              out_file <<rho(j,i)<<"\n";
          }
      }
        
      out_file << "VECTORS "<<"U "<<"double\n";

      
      for(int j=0; j < Ny; j++)
      {
          for(int i=0; i <Nx; i++)
          {
              out_file <<ux(j,i)<<" "<<uy(j,i)<<" 0\n";
          }
      }

      out_file.close();
 
    }
    else
    {
        std::cout<<"Something went wrong on line 260 in Solver.\n";
    }
  }

  void output_VTK(int s,int plot_every)
  {
    if(std::is_same_v< DeviceType, TNL::Devices::Cuda > )
    {
      //std::cout<<"\nCuda -> Host output\n"<<std::endl;
      ArrayType2D_Host rho_out;
      ArrayType2D_Host ux_out;
      ArrayType2D_Host uy_out;
      
      rho_out = rho;
      ux_out = ux;
      uy_out = uy; 
      std::string step =std::to_string(s/plot_every);

      std::ofstream out_file("results/LBE."+step+".vtk");

      out_file << "# vtk DataFile Version 2.0\n" ;
      out_file << "LBE mesh\n" ;
      out_file << "ASCII\n";
      out_file << "DATASET STRUCTURED_POINTS\n" ;
      out_file << "DIMENSIONS "<<Nx<<" "<<Ny<<" 1\n" ;
      out_file << "ASPECT_RATIO 1 1 1\n";
      out_file << "ORIGIN 0 0 0\n";
      out_file << "POINT_DATA "<<Nx*Ny<<"\n";

      out_file << "SCALARS "<<"rho "<<"double 1\n";
      out_file << "LOOKUP_TABLE default\n";

      RealType Crho = Cm/(Cl*Cl*Cl);
      for(int j=0; j < Ny; j++)
      {
          for(int i=0; i < Nx; i++)
          {
              out_file <<rho_out(j,i)*Crho<<"\n";
          }
      }
        
      out_file << "VECTORS "<<"U "<<"double\n";

      
      for(int j=0; j < Ny; j++)
      {
          for(int i=0; i <Nx; i++)
          {
              out_file <<ux_out(j,i)*Cu<<" "<<uy_out(j,i)*Cu<<" 0\n";
          }
      }

      out_file.close();
    }

    else if (std::is_same<DeviceType, TNL::Devices::Host >::value )
    {
      //std::cout<<"\nHost output\n"<<std::endl;
      std::string step =std::to_string(s/plot_every);

      std::ofstream out_file("results/LBE."+step+".vtk");

      out_file << "# vtk DataFile Version 2.0\n" ;
      out_file << "LBE mesh\n" ;
      out_file << "ASCII\n";
      out_file << "DATASET STRUCTURED_POINTS\n" ;
      out_file << "DIMENSIONS "<<Nx<<" "<<Ny<<" 1\n" ;
      out_file << "ASPECT_RATIO 1 1 1\n";
      out_file << "ORIGIN 0 0 0\n";
      out_file << "POINT_DATA "<<Nx*Ny<<"\n";

      out_file << "SCALARS "<<"rho "<<"double 1\n";
      out_file << "LOOKUP_TABLE default\n";

      RealType Crho = Cm/(Cl*Cl*Cl);
      for(int j=0; j < Ny; j++)
      {
          for(int i=0; i < Nx; i++)
          {
              out_file <<rho(j,i)*Crho<<"\n";
          }
      }
        
      out_file << "VECTORS "<<"U "<<"double\n";

      
      for(int j=0; j < Ny; j++)
      {
          for(int i=0; i <Nx; i++)
          {
              out_file <<ux(j,i)*Cu<<" "<<uy(j,i)*Cu<<" 0\n";
          }
      }

      out_file.close();
    }

  }

  void convert_to_lattice(RealType L_fyz, RealType U_fyz, RealType rho_fyz, RealType ny_fyz, RealType U_lb)
  {
    // prenasobim lattice jednotky a dostanu fyikalni
    // vydelim fyz a dostanu lattice
    if(U_lb > 0.1)
    {
        std::cout<<"\n Lattice speed should not be higher than 0.1 due to stability close to speed of sound.\n";        
    }

    assert(U_lb<cs);

    Cl = L_fyz/Nx;
    Cu = U_fyz/U_lb;
    Cu_inverse = 1/Cu;
    Ct = Cl/Cu;
    Ct_pub = Ct;
    RealType Crho = rho_fyz;
    Cm = Crho*Cl*Cl*Cl;

    RealType Re = U_fyz*L_fyz/ny_fyz;

    std::cout<<"\n Re is "<<Re<<"\n";

    ny = ny_fyz*Ct/Cl/Cl;

    tau = ny/3 + 0.5f;

    if(tau < 0.51)
    {
        std::cout<<"Tau is too small. Consider higher resolution or higher lattice speed. \n";
        std::cout<<"Tau is "<<tau<<"\n";
    }
    
    else if(tau > 0.99 )
    {
        std::cout<<"Tau is too high. Consider lower resolution or lower lattice speed. \n";
        std::cout<<"Tau is "<<tau<<"\n";
    }

    else
    {
        std::cout<<"\n Tau is "<<tau<<" <-(0.5;1)\n";
    
    }

    assert(tau>0.5 && tau<0.99);

    omega_0=1/tau;

    std::cout<<"\n Conversion undergone successfully."<<"\n";

  }

  ArrayType2D rho;
  ArrayType2D ux;
  ArrayType2D uy;
  ArrayType2D omega;

  RealType Ct_pub;
  RealType err = 1.f;

private:
  //Model parameters
  /*
  6    2    5
    \  |  /
  3 -- 0 -- 1 
    /  |  \
  7    4    8
  */

  // view df = this-> get.view

  const int Nvel = 9 ; // Q number of possible descrete velocities
  const int cx_pos[9] = { 0, 1, 0, -1, 0, 1, -1, -1,  1};
  const int cy_pos[9] = { 0, 0, 1, 0, -1, 1,  1, -1, -1};
  const int c[9][2]={
          {0,0},
          {1,0},
          {0,1},
          {-1,0},
          {0,-1},
          {1,1},
          {-1,1},
          {-1,-1},
          {1,-1}};

  const int c_rever[9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };
  const RealType weight[9] ={4.f/9.f,1.f/9.f,1.f/9.f,1.f/9.f,1.f/9.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f};
  const RealType cs = 1/sqrt(3.f);
  const RealType cs4 = 1/9.f;
  const RealType CLES = 1.73f;

  RealType ny;
  RealType tau;
  RealType omega_0;

  
  ArrayType3D df;
  ArrayType3D df_post;

  ArrayType2D mesh;
  ArrayType2D velocities_x;
  ArrayType2D velocities_y;
  ArrayType2D ux0er;
  ArrayType2D uy0er;

  int Nx;
  int Ny;
  int s;

  RealType Fx;
  RealType Fy;
  RealType g;

  RealType Cl;
  RealType Ct;    
  RealType Cm;
  RealType Cu_inverse;
  RealType Cu;


};

#endif