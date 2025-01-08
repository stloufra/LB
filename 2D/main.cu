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

using DeviceType = TNL::Devices::Cuda;
using DeviceTypeHost = TNL::Devices::Host;

using RealType = float;

int main()
{
    const RealType L = 0.1f;                //[m] - x dimension
    const int Nx = 1000;                    //[1]
    const int Ny = 150;                     //[1]

    const RealType rho=1000.f;              //[kg/m3]
    const RealType ny=1e-5f;                //[m2/s]

    const RealType ux=0.01f;                //[m/s] // 0.01 ok 0.1 fail
    const RealType ux_guess=0.02f;          //[m/s]
    const RealType uy=0.f;                  //[m/s]
    const RealType u_max_lattice =0.09f;    //[0]

    const RealType Fx = 0.0f;               //[kg/m2/s2]  <- force density (3rd dimension in 2D is equal to 1)
    const RealType Fy = 0.0f;               //[kg/m2/s2]  <- force density (3rd dimension in 2D is equal to 1)

    const RealType time =5.f;              //[s]

    const int plot_every = 10.f;            //[s]
    const int err_every_it = 1000;          //[it]

    int iterations;
    
    Mesher<RealType, DeviceTypeHost > mesh_rectangle(Ny,Nx);   

    //objects

    Obj_rectangle lower_wall( -1, Nx , -1, -1);
    Obj_rectangle symmetry( 0, Nx-1 , 0, 0);
    Obj_rectangle upper_wall( -1, Nx, Ny , Ny );
    Obj_rectangle front_wall(-1, -1, 0, Ny-1 );
    Obj_rectangle back_wall(Nx , Nx, 0, Ny-1);
    Obj_rectangle inlet(0, 0, 0, Ny-1 );
    Obj_rectangle outlet(Nx-1 , Nx-1, 0, Ny-1);
    Obj_cylinder cylinder(Ny/5, Nx/4, Ny/2);

    // MESH - structured bolean values of BC
    // 0 = solid | 1 = fluid | 2 = equilibrium inlet (ux, left) | 3 = outlet (rho=1, right) | 4 = moving wall (up) | 5 = moving wall (down) | 6 = symmetry down

    mesh_rectangle.meshing(lower_wall,0);
    mesh_rectangle.meshing(upper_wall,0);
    mesh_rectangle.meshing(front_wall,0);
    mesh_rectangle.meshing(back_wall,0);

    mesh_rectangle.meshing(symmetry,6);
    mesh_rectangle.meshing(outlet, 3);
    mesh_rectangle.meshing_moving(inlet, ux, 0, 2);

    mesh_rectangle.meshing(cylinder, 0);


    //output mesh
    mesh_rectangle.output_VTK();

    Solver_sac<RealType, DeviceType> solver(Ny,Nx,mesh_rectangle);

    //non-dimensionalize
    solver.convert_to_lattice(L, ux_guess, rho, ny, u_max_lattice);

    int plot_every_it = std::ceil(plot_every/solver.Ct_pub);
    std::cout<<"\nPlotting every " << plot_every_it << " iterations.\n";
    iterations = std::ceil(time/solver.Ct_pub);
    std::cout<<"\nCalculation will run for "<<iterations<<" iterations.\n";

    solver.initialization_eq(rho, ux, 0.0001f, Fx, Fy,0);

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
    while(k<iterations) //err>=10e-4)
    {
        k++;

        
        timer_collision.start();
        solver.collision_SRT();
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
            solver.appendError(k);
            printf("\n err=%e, k=%d \n" ,solver.err, k);


            if (std::isnan(solver.err)) {
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
    logger.writeSystemInformation(true);
    logger.writeHeader("Loop");
    timer_loop.writeLog( logger, 0 );
    logger.writeSeparator();
    logger.writeHeader("Collision");
    timer_collision.writeLog( logger, 0 );
    logger.writeSeparator();
    logger.writeHeader("Streaming");
    timer_streaming.writeLog( logger, 0 );
    logger.writeSeparator();
    logger.writeHeader("Boundary Conditions");
    timer_bounceback.writeLog( logger, 0 );
    logger.writeSeparator();
    logger.writeHeader("Moments");
    timer_postpro.writeLog( logger, 0 );
    logger.writeSeparator();
    logger.writeHeader("Error Calculation");
    timer_err.writeLog( logger, 0 );
    logger.writeSeparator();
    logger.writeHeader("Output");
    timer_output.writeLog( logger, 0 );


    return 0;
}

