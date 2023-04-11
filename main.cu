#include <iostream>
#include <fstream>
#include "Mesher.h"
#include "Solver_sac.h"
#include "Obj_cylinder.h"
#include "Obj_rectangle.h"
#include "Obj_template.h"
#include <TNL/Timer.h>
#include <TNL/Logger.h>

using DeviceType = TNL::Devices::Host;
using DeviceTypeHost = TNL::Devices::Host;

using RealType = double;



int main()
{
    const double L = 0.1;               //[m]
    const int Nx = 300;                 //[1]
    const int Ny = 100;                 //[1]

    const double rho=1000;              //[kg/m3]
    const double ny=10e-5;              //[m2/s]

    const double ux=0.100;              //[m/s]
    const double ux_guess=0.1;          //[m/s]
    const double uy=0.000;              //[m/s]
    const double u_max_lattice =0.1;  //[0]

    const double Fx = 10;               //[kg/m2/s2]  <- force density (3rd dimension in 2D is equal to 1)
    const double Fy = 0.0;              //[kg/m2/s2]  <- force density (3rd dimension in 2D is equal to 1)

    const double time =1;               //[s]  
    const double plot_every=0.01;       //[s]

    int plot_every_it;
    int iterations;

    
    
    Mesher<RealType, DeviceTypeHost > mesh_rectangle(Ny,Nx);   

    //objects !! pristup od 0 !! horní index o 1 mensi je to 

    Obj_rectangle lower_wall( -1.0, Nx , -1.0, -1.0);
    Obj_rectangle upper_wall( -1.0, Nx, Ny , Ny );
    Obj_rectangle inlet(-1, -1, 0, Ny-1 );
    Obj_rectangle outlet(Nx , Nx, 0, Ny-1);
    Obj_cylinder cylinder(Ny/5, Nx/4,Ny/2+0.05*Ny);

    // MESH - structured bolean values of BC
    // 0 = solid | 1 = fluid | 2 = primitive inlet vertical | 3 = outlet (rho=1, right) | 4 = moving wall up | 5 = moving wall down | 6 = outlet (rh=1, left)

    mesh_rectangle.meshing(lower_wall,0);
    mesh_rectangle.meshing(upper_wall,0);
    mesh_rectangle.meshing(cylinder, 0);
    mesh_rectangle.meshing(outlet, 3);
    mesh_rectangle.meshing_moving(inlet, ux, 0, 2);

    //output mesh
    mesh_rectangle.output_VTK();

    Solver_sac<RealType, DeviceType> solver(Ny,Nx,mesh_rectangle);
    solver.convert_to_lattice(L, ux_guess, rho, ny, u_max_lattice);

    
    plot_every_it = std::ceil(plot_every/solver.Ct_pub);
    std::cout<<"\nPlotting every " << plot_every_it << " iterations.\n";
    iterations = std::ceil(time/solver.Ct_pub);
    std::cout<<"\nCalculation will run for "<<iterations<<" iterations.\n";
    
    solver.initialization_eq(rho, ux, uy, Fx, Fy, 0);

    solver.output_VTK_lattice();
    solver.output_VTK(0,plot_every_it);
    
    //solver run
    
    Timer timer;
    Logger logger(50, std::cout);

    timer.start();

    int k = 0;
    while(k<iterations) //err>=10e-4)
    {
        k++;
        solver.collision();
        solver.streaming();
        solver.bounce_back();
        solver.postpro();

        if(k%500==0 && k!=0)
        {
            //solver.Err();
            //printf("\n err=%e ux_center=%e uy_center=%e rho_center=%e k=%d\n",solver.err,solver.ux(Ny/2,Nx/2),solver.uy(Ny/2,Nx/2),solver.rho(Ny/2,Nx/2), k);
            solver.output_VTK_lattice();
             if (std::isnan(solver.err)) 
             {
                std::cout << "\n Error is NaN, breaking out.\n";
                break;
             }
        }

        

        if(k%plot_every_it==0)
        {
            solver.output_VTK(k,plot_every_it);
        }
   
    }

    timer.stop();
    
    timer.writeLog( logger, 0 );

    return 0;
}

