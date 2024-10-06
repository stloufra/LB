#include <iostream>
#include <fstream>
#include "./src/geo/Mesher.h"
#include "./src/sol/Solver_sac.h"
#include "./src/sol/Solver_LES.h"
#include "./src/geo/Obj_cylinder.h"
#include "./src/geo/Obj_rectangle.h"
#include "./src/geo/Obj_template.h"
#include <TNL/Timer.h>
#include <TNL/Logger.h>

using DeviceType = TNL::Devices::Cuda; //sac
using DeviceTypeHost = TNL::Devices::Host;

using RealType = float;



int main()
{
    const RealType L = 15.0f;                //[m]
    const int Nx = 3000;                     //[1]
    const int Ny = 400;                     //[1]

    const RealType rho=1000.f;              //[kg/m3]
    const RealType ny=1e-6f;               //[m2/s]

    const RealType ux=1e-4f;               //[m/s] // 0.01 ok 0.1 fail
    const RealType ux_guess=2e-4f;          //[m/s]
    const RealType uy=0.f;                  //[m/s]
    const RealType u_max_lattice =0.09f;    //[0]

    const RealType Fx = 0.0f;               //[kg/m2/s2]  <- force density (3rd dimension in 2D is equal to 1)
    const RealType Fy = 0.0f;               //[kg/m2/s2]  <- force density (3rd dimension in 2D is equal to 1)

    const RealType time =10000000.f;               //[s]
    const RealType plot_every=200.f;         //[s]

    int plot_every_it;
    int iterations;

    bool run_simulation = true;
    
    
    Mesher<RealType, DeviceTypeHost > mesh_rectangle(Ny,Nx);   

    //objects !! pristup od 0 !! horní index o 1 mensi je to 

    Obj_rectangle lower_wall( -1.f, Nx , -1.f, -1.f);
    Obj_rectangle upper_wall( -1.f, Nx, Ny , Ny );
    Obj_rectangle front_wall(-1, -1, 0, Ny-1 );
    Obj_rectangle back_wall(Nx , Nx, 0, Ny-1);
    Obj_rectangle inlet(0, 0, 0, Ny-1 );
    Obj_rectangle outlet(Nx-1 , Nx-1, 0, Ny-1);
    //Obj_cylinder cylinder(Ny/5, Nx/4,Ny/2+0.05f*Ny);
    Obj_rectangle blockage(-1, Nx/5, -1, 3*Ny/5);

    // MESH - structured bolean values of BC
    // 0 = solid | 1 = fluid | 2 = primitive inlet vertical | 3 = outlet (rho=1, right) | 4 = moving wall up | 5 = moving wall down | 6 = outlet (rh=1, left)

    mesh_rectangle.meshing(lower_wall,0);
    mesh_rectangle.meshing(upper_wall,0);
    mesh_rectangle.meshing(front_wall,0);
    mesh_rectangle.meshing(back_wall,0);
    mesh_rectangle.meshing(outlet, 3);
    mesh_rectangle.meshing_moving(inlet, ux, 0, 2);

    //mesh_rectangle.meshing(cylinder, 0);
    mesh_rectangle.meshing(blockage,0);

    //output mesh
    mesh_rectangle.output_VTK();

    Solver_sac<RealType, DeviceType> solver(Ny,Nx,mesh_rectangle);
    solver.convert_to_lattice(L, ux_guess, rho, ny, u_max_lattice);

    
    plot_every_it = std::ceil(plot_every/solver.Ct_pub);
    std::cout<<"\nPlotting every " << plot_every_it << " iterations.\n";
    iterations = std::ceil(time/solver.Ct_pub);
    std::cout<<"\nCalculation will run for "<<iterations<<" iterations.\n";

    plot_every_it = 500;
    
    solver.initialization_eq(rho, ux, uy, Fx, Fy, 0);

    solver.output_VTK_lattice();
    solver.output_VTK(0,plot_every_it);


    //solver run
    
    Timer timer_loop;
    Timer timer_collision;
    Timer timer_streaming;
    Timer timer_bounceback;
    Timer timer_postpro;
    Timer timer_err;
    Timer timer_output;

    Logger logger(50, std::cout);

    timer_loop.start();


    int k = 0;
    while(k<iterations && run_simulation) //err>=10e-4)
    {
        k++;
        /*timer_postpro.start();
        solver.omegaLES();
        timer_postpro.stop();*/
        
        timer_collision.start();
        solver.collision();
        timer_collision.stop();


        timer_streaming.start();
        solver.streaming();
        timer_streaming.stop();

        timer_bounceback.start();
        solver.bounce_back();
        timer_bounceback.stop();

        timer_postpro.start();
        solver.postpro();
        timer_postpro.stop();

        if(k%500==0 && k!=0)
        {

            timer_err.start();
            solver.Err();
            //printf("\n err=%e ux_center=%e uy_center=%e rho_center=%e k=%d\n",solver.err,solver.ux.getView()(Ny/2,Nx/2),solver.uy.getView()(Ny/2,Nx/2),solver.rho.getView()(Ny/2,Nx/2), k);
            printf("\n err=%e, k=%d \n" ,solver.err, k);
            if (std::isnan(solver.err))
             {
                std::cout << "\n Error is NaN, breaking out.\n";
                break;
             }

            timer_err.stop();
        }

        

        if(k%plot_every_it==0)
        {   
            timer_output.start();
            solver.output_VTK(k,plot_every_it);
            timer_output.stop();
        }
   
    }

    timer_loop.stop();
    
    logger.writeHeader("Timing of sections");
    //logger.writeSystemInformation(true);
    timer_loop.writeLog( logger, 0 );
    logger.writeSeparator();
    timer_collision.writeLog( logger, 0 );
    logger.writeSeparator();
    timer_streaming.writeLog( logger, 0 );
    logger.writeSeparator();
    timer_bounceback.writeLog( logger, 0 );
    logger.writeSeparator();
    timer_postpro.writeLog( logger, 0 );
    logger.writeSeparator();
    timer_err.writeLog( logger, 0 );
    logger.writeSeparator();
    timer_output.writeLog( logger, 0 );


    return 0;
}

