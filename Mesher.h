#ifndef MESHER_H
#define MESHER_H

#pragma once

#include <fstream>
#include <iostream>
#include <omp.h>
#include <cassert>
#include "Obj_template.h"
#include <TNL/Containers/NDArray.h>


using namespace TNL;
using namespace TNL::Containers;

template< typename RealType, typename DeviceType >
class Mesher
{

using ArrayType =  TNL::Containers::NDArray< RealType,
                                             SizesHolder<int, 0, 0 >,
                                             std::index_sequence< 0, 1 > ,
                                             DeviceType>;

public:
    Mesher()=delete;


    Mesher(int Ny_in, int Nx_in):
        mesh(),
        velocities_x(),
        velocities_y()
        
    {   
        Nx = Nx_in;
        Ny = Ny_in;

        mesh.setSizes(Ny+2, Nx+2);
        velocities_x.setSizes(Ny+2, Nx+2);
        velocities_y.setSizes(Ny+2, Nx+2);
        

        for(int j=-1; j < Ny + 1; j++)
        {
            for(int i=-1; i < Nx +1 ; i++)
            {
                mesh(j+1,i+1) = 1.0;
                velocities_x(j+1,i+1) = 0.0;
                velocities_y(j+1,i+1) = 0.0;
            }
        }

    }

    ~Mesher() {}


//mesh functions
    void meshing(Obj_template &obj, int type)
    {
        for(int j=-1; j < Ny + 1; j++)
        {
            for(int i=-1; i < Nx +1 ; i++)
            {
                if(obj.is_inside(i,j))
                {
                    mesh(j+1,i+1) = type;
                }
            }  
        }

        std::cout<< "Meshing of ";
        std::cout<< PATCH(type);
        std::cout <<" undergone successfully."<<"\n";
    }

    void meshing_moving(Obj_template &obj, double v_x, double v_y, int type)
    {
        for(int j=-1; j < Ny + 1; j++)
        {
            for(int i=-1; i < Nx + 1 ; i++)
            { 
                if(obj.is_inside(i,j))
                {
                    mesh(j+1,i+1) = type;
                    velocities_x(j+1,i+1)= v_x;
                    velocities_y(j+1,i+1) = v_y;
                }
            }  
        }

        std::cout<< "Meshing of ";
        std::cout<< PATCH(type);
        std::cout <<" undergone successfully."<<"\n";
    }

    //help functions
    void output_VTK()
    {

        int Nx2 =  Nx+2 ;
        int Ny2 =  Ny+2 ;

        std::ofstream out_file("results/Mesh.vtk");

        out_file << "# vtk DataFile Version 2.0\n" ;
        out_file << "LBE mesh\n" ;
        out_file << "ASCII\n";
        out_file << "DATASET STRUCTURED_POINTS\n" ;
        out_file << "DIMENSIONS "<<Nx2<<" "<<Ny2<<" 1\n" ;
        out_file << "ASPECT_RATIO 1 1 1\n";
        out_file << "ORIGIN 0 0 0\n";
        out_file << "POINT_DATA "<<Nx2*Ny2<<"\n";

        out_file << "SCALARS "<<"velocities_x "<<"double 1\n";
        out_file << "LOOKUP_TABLE default\n";
        for(int j=-1; j < Ny + 1; j++)
        {
            for(int i=-1; i < Nx +1 ; i++)
            { 
            out_file << velocities_x(j+1,i+1) << "\n";
            }
        }


        out_file << "SCALARS "<<"velocities_y "<<"double 1\n";
        out_file << "LOOKUP_TABLE default\n";
        for(int j=-1; j < Ny + 1; j++)
        {
            for(int i=-1; i < Nx +1 ; i++)
            {
            out_file <<velocities_y(j+1,i+1)<<"\n";
            }
        }

        out_file << "SCALARS "<<"mesh "<<"int 1\n";
        out_file << "LOOKUP_TABLE default\n";
        for(int j=-1; j < Ny + 1; j++)
        {
            for(int i=-1; i < Nx +1 ; i++)
            {
            out_file << mesh(j+1,i+1) <<"\n";
            }
        }
    }
    


    ArrayType mesh;
    ArrayType velocities_x;
    ArrayType velocities_y;
    
private:

    int Nx;
    int Ny;

    enum PATCH 
    {
        solid = 0,
        fluid = 1,
        inlet_right_vertical = 2,
        outlet_equlibrium = 3,
        wall_moving_up = 4,
        wall_moving_down = 5
    }; 
};

#endif