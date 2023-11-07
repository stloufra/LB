//
// Created by stloufra on 10/16/23.
//

#include <fstream>
#include <iostream>
#include <string>
#include "../traits/LBMTraits.h"

#ifndef INC_24_4_GEO_TO_TNL_OUTPUTER_H
#define INC_24_4_GEO_TO_TNL_OUTPUTER_H

class outputerVTK
{
    using RealType = LBMTraits::RealType;
    using DeviceType = LBMTraits::DeviceType;
    using VectorType = LBMTraits::VectorType;
    using LBMDataPointer = TNL::Pointers::SharedPointer<LBMData, DeviceType>;
    using LBMConstantsPointer = TNL::Pointers::SharedPointer<LBMConstants, DeviceType>;
    using ArrayTypeVariablesScalarHost = LBMTraits::ArrayTypeVariablesScalarHost;
    using ArrayTypeVariablesVectorHost = LBMTraits::ArrayTypeVariablesVectorHost;
    using ArrayTypeDFunctionHost = LBMTraits::ArrayTypeDFunctionHost;
public:

    static void MeshVTK(LBMDataPointer &Data, LBMConstantsPointer &Constants, std::string filename) {
        std::cout << "\nWriting vtk file - " << filename << ".vtk\n";

        std::ofstream out_file( "results/" + filename + ".vtk");

        out_file << "# vtk DataFile Version 2.0\n";
        out_file << "LBE " << filename << "\n";
        out_file << "ASCII\n";
        out_file << "DATASET STRUCTURED_POINTS\n";
        out_file << "DIMENSIONS " << Constants->dimX_int << " " << Constants->dimY_int << " " << Constants->dimZ_int
                 << "\n";
        out_file << "ASPECT_RATIO " << 1.f / Constants->resolution_factor << " " << 1.f / Constants->resolution_factor
                 << " "
                 << 1.f / Constants->resolution_factor << "\n";
        out_file << "ORIGIN " << Constants->BBminx - Constants->additional_factor << " "
                 << Constants->BBminy - Constants->additional_factor << " "
                 << Constants->BBminz - Constants->additional_factor << "\n";
        out_file << "POINT_DATA " << Constants->dimX_int * Constants->dimY_int * Constants->dimZ_int << "\n";

        out_file << "SCALARS " << "object " << "double 1\n";
        out_file << "LOOKUP_TABLE default\n";

        for (int k = 0; k < Constants->dimZ_int; k++) {
            for (int j = 0; j < Constants->dimY_int; j++) {
                for (int i = 0; i < Constants->dimX_int; i++) {
                    out_file << Data->meshFluidHost(i, j, k) << "\n";
                }
            }
        }

        std::cout << "Mesh written.\n";
    }

    static void variablesLatticeVTK(LBMDataPointer &Data, LBMConstantsPointer &Constants, int step, bool verbose) {

        if (std::is_same_v<DeviceType, TNL::Devices::Cuda>) {
            if (verbose) { std::cout << "\nCuda -> Host output Lattice units\n" << std::endl; }

            ArrayTypeVariablesScalarHost rho_out;
            ArrayTypeVariablesScalarHost p_out;
            ArrayTypeVariablesVectorHost u_out;


            rho_out = Data->rho;
            u_out = Data->u;
            p_out = Data->p;


            std::ofstream out_file("results/variablesLattice"+std::to_string(step)+".vtk");


            out_file << "# vtk DataFile Version 2.0\n";
            out_file << "LBE variables\n";
            out_file << "ASCII\n";
            out_file << "DATASET STRUCTURED_POINTS\n";
            out_file << "DIMENSIONS " << Constants->dimX_int << " " << Constants->dimY_int << " " << Constants->dimZ_int
                     << "\n";
            out_file << "ASPECT_RATIO " << 1.f / Constants->resolution_factor << " "
                     << 1.f / Constants->resolution_factor
                     << " "
                     << 1.f / Constants->resolution_factor << "\n";
            out_file << "ORIGIN " << Constants->BBminx - Constants->additional_factor << " "
                     << Constants->BBminy - Constants->additional_factor << " "
                     << Constants->BBminz - Constants->additional_factor << "\n";
            out_file << "POINT_DATA " << Constants->dimX_int * Constants->dimY_int * Constants->dimZ_int << "\n";

            out_file << "SCALARS " << "rho " << "double 1\n";
            out_file << "LOOKUP_TABLE default\n";

            for (int k = 0; k < Constants->dimZ_int; k++) {
                for (int j = 0; j < Constants->dimY_int; j++) {
                    for (int i = 0; i < Constants->dimX_int; i++) {
                        out_file << rho_out(i, j, k) << "\n";
                    }
                }
            }

            out_file << "SCALARS " << "pressure " << "double 1\n";
            out_file << "LOOKUP_TABLE default\n";

            for (int k = 0; k < Constants->dimZ_int; k++) {
                for (int j = 0; j < Constants->dimY_int; j++) {
                    for (int i = 0; i < Constants->dimX_int; i++) {
                        out_file << p_out(i, j, k) << "\n";
                    }
                }
            }


            out_file << "VECTORS " << "U " << "double\n";

            for (int k = 0; k < Constants->dimZ_int; k++) {
                for (int j = 0; j < Constants->dimY_int; j++) {
                    for (int i = 0; i < Constants->dimX_int; i++) {
                        out_file << u_out(i, j, k).x() << " " << u_out(i, j, k).y() << " "
                                 << u_out(i, j, k).z() << "\n";
                    }
                }
            }


            out_file.close();

        } else if (std::is_same<DeviceType, TNL::Devices::Host>::value) {
            if (verbose) { std::cout << "\nHost -> Host output Lattice units\n" << std::endl; }

            std::ofstream out_file("results/variablesLattice"+std::to_string(step)+".vtk");


            out_file << "# vtk DataFile Version 2.0\n";
            out_file << "LBE variables\n";
            out_file << "ASCII\n";
            out_file << "DATASET STRUCTURED_POINTS\n";
            out_file << "DIMENSIONS " << Constants->dimX_int << " " << Constants->dimY_int << " " << Constants->dimZ_int
                     << "\n";
            out_file << "ASPECT_RATIO " << 1.f / Constants->resolution_factor << " "
                     << 1.f / Constants->resolution_factor
                     << " "
                     << 1.f / Constants->resolution_factor << "\n";
            out_file << "ORIGIN " << Constants->BBminx - Constants->additional_factor << " "
                     << Constants->BBminy - Constants->additional_factor << " "
                     << Constants->BBminz - Constants->additional_factor << "\n";
            out_file << "POINT_DATA " << Constants->dimX_int * Constants->dimY_int * Constants->dimZ_int << "\n";

            out_file << "SCALARS " << "rho " << "double 1\n";
            out_file << "LOOKUP_TABLE default\n";

            for (int k = 0; k < Constants->dimZ_int; k++) {
                for (int j = 0; j < Constants->dimY_int; j++) {
                    for (int i = 0; i < Constants->dimX_int; i++) {
                         out_file << Data->rho(i, j, k) << "\n";
                    }
                }
            }

            out_file << "SCALARS " << "pressure " << "double 1\n";
            out_file << "LOOKUP_TABLE default\n";

            for (int k = 0; k < Constants->dimZ_int; k++) {
                for (int j = 0; j < Constants->dimY_int; j++) {
                    for (int i = 0; i < Constants->dimX_int; i++) {
                        out_file << Data->p(i, j, k) << "\n";
                    }
                }
            }

            out_file << "VECTORS " << "U " << "double\n";

            for (int k = 0; k < Constants->dimZ_int; k++) {
                for (int j = 0; j < Constants->dimY_int; j++) {
                    for (int i = 0; i < Constants->dimX_int; i++) {
                        out_file << Data->u(i, j, k).x() << " " << Data->u(i, j, k).y() << " " << Data->u(i, j, k).z()
                                 << "\n";
                    }
                }
            }

            out_file.close();

        }
    }

    static void distributionFunctionVTK(LBMDataPointer &Data, LBMConstantsPointer &Constants, int velocity, bool verbose) {

        if (std::is_same_v<DeviceType, TNL::Devices::Cuda>) {
            if (verbose) { std::cout << "\nCuda -> Host output distribution function\n" << std::endl; }


            ArrayTypeDFunctionHost df_out;


            df_out = Data->df;


            std::ofstream out_file("results/distributionFunction"+ std::to_string(velocity)+".vtk");


            out_file << "# vtk DataFile Version 2.0\n";
            out_file << "LBE variables\n";
            out_file << "ASCII\n";
            out_file << "DATASET STRUCTURED_POINTS\n";
            out_file << "DIMENSIONS " << Constants->dimX_int << " " << Constants->dimY_int << " " << Constants->dimZ_int
                     << "\n";
            out_file << "ASPECT_RATIO " << 1.f / Constants->resolution_factor << " "
                     << 1.f / Constants->resolution_factor
                     << " "
                     << 1.f / Constants->resolution_factor << "\n";
            out_file << "ORIGIN " << Constants->BBminx - Constants->additional_factor << " "
                     << Constants->BBminy - Constants->additional_factor << " "
                     << Constants->BBminz - Constants->additional_factor << "\n";
            out_file << "POINT_DATA " << Constants->dimX_int * Constants->dimY_int * Constants->dimZ_int << "\n";

            out_file << "SCALARS " << "DFunciton_" + std::to_string(velocity) << " double 1\n";
            out_file << "LOOKUP_TABLE default\n";

            for (int k = 0; k < Constants->dimZ_int; k++) {
                for (int j = 0; j < Constants->dimY_int; j++) {
                    for (int i = 0; i < Constants->dimX_int; i++) {
                        out_file << df_out(i, j, k, velocity) << "\n";
                    }
                }
            }


            out_file.close();

        } else if (std::is_same<DeviceType, TNL::Devices::Host>::value) {
            if (verbose) { std::cout << "\nHost -> Host output distribution function.\n" << std::endl; }

            std::ofstream out_file("results/distributionFunction"+ std::to_string(velocity)+".vtk");


            out_file << "# vtk DataFile Version 2.0\n";
            out_file << "LBE variables\n";
            out_file << "ASCII\n";
            out_file << "DATASET STRUCTURED_POINTS\n";
            out_file << "DIMENSIONS " << Constants->dimX_int << " " << Constants->dimY_int << " " << Constants->dimZ_int
                     << "\n";
            out_file << "ASPECT_RATIO " << 1.f / Constants->resolution_factor << " "
                     << 1.f / Constants->resolution_factor
                     << " "
                     << 1.f / Constants->resolution_factor << "\n";
            out_file << "ORIGIN " << Constants->BBminx - Constants->additional_factor << " "
                     << Constants->BBminy - Constants->additional_factor << " "
                     << Constants->BBminz - Constants->additional_factor << "\n";
            out_file << "POINT_DATA " << Constants->dimX_int * Constants->dimY_int * Constants->dimZ_int << "\n";

            out_file << "SCALARS " << "DFunciton_" + std::to_string(velocity) << " double 1\n";
            out_file << "LOOKUP_TABLE default\n";

            for (int k = 0; k < Constants->dimZ_int; k++) {
                for (int j = 0; j < Constants->dimY_int; j++) {
                    for (int i = 0; i < Constants->dimX_int; i++) {
                        out_file << Data -> df(i, j, k, velocity) << "\n";
                    }
                }
            }


            out_file.close();

        }
    }

    static void variablesVTK(LBMDataPointer &Data, LBMConstantsPointer &Constants, int s, int plot_every, bool verbose) {

        if (std::is_same_v<DeviceType, TNL::Devices::Cuda>) {
            if (verbose) { std::cout << "\nCuda -> Host output Lattice units\n" << std::endl; }



            Data->rho_out = Data->rho;
            Data->u_out = Data->u;
            Data->p_out = Data->p;

            std::string step =std::to_string(s/plot_every);

            std::ofstream out_file("results/variables."+step+".vtk");


            out_file << "# vtk DataFile Version 2.0\n";
            out_file << "LBE variables\n";
            out_file << "ASCII\n";
            out_file << "DATASET STRUCTURED_POINTS\n";
            out_file << "DIMENSIONS " << Constants->dimX_int << " " << Constants->dimY_int << " " << Constants->dimZ_int
                     << "\n";
            out_file << "ASPECT_RATIO " << 1.f / Constants->resolution_factor << " "
                     << 1.f / Constants->resolution_factor
                     << " "
                     << 1.f / Constants->resolution_factor << "\n";
            out_file << "ORIGIN " << Constants->BBminx - Constants->additional_factor << " "
                     << Constants->BBminy - Constants->additional_factor << " "
                     << Constants->BBminz - Constants->additional_factor << "\n";
            out_file << "POINT_DATA " << Constants->dimX_int * Constants->dimY_int * Constants->dimZ_int << "\n";

            out_file << "SCALARS " << "rho " << "double 1\n";
            out_file << "LOOKUP_TABLE default\n";

            for (int k = 0; k < Constants->dimZ_int; k++) {
                for (int j = 0; j < Constants->dimY_int; j++) {
                    for (int i = 0; i < Constants->dimX_int; i++) {
                        out_file << Data->rho_out(i, j, k) * Constants->Crho << "\n";
                    }
                }
            }

            out_file << "SCALARS " << "p " << "double 1\n";
            out_file << "LOOKUP_TABLE default\n";

            for (int k = 0; k < Constants->dimZ_int; k++) {
                for (int j = 0; j < Constants->dimY_int; j++) {
                    for (int i = 0; i < Constants->dimX_int; i++) {
                        out_file << Data->p_out(i, j, k) * Constants->Cpressure << "\n";
                    }
                }
            }


            out_file << "VECTORS " << "U " << "double\n";

            for (int k = 0; k < Constants->dimZ_int; k++) {
                for (int j = 0; j < Constants->dimY_int; j++) {
                    for (int i = 0; i < Constants->dimX_int; i++) {
                        out_file << Data->u_out(i, j, k).x() * Constants->Cu << " "
                                 << Data->u_out(i, j, k).y() * Constants->Cu << " "
                                 << Data->u_out(i, j, k).z() * Constants->Cu << "\n";
                    }
                }
            }


            out_file.close();

        } else if (std::is_same<DeviceType, TNL::Devices::Host>::value) {
            if (verbose) { std::cout << "\nHost output Lattice units\n" << std::endl; }

            std::string step =std::to_string(s/plot_every);

            std::ofstream out_file("results/variables."+step+".vtk");

            out_file << "# vtk DataFile Version 2.0\n";
            out_file << "LBE variables\n";
            out_file << "ASCII\n";
            out_file << "DATASET STRUCTURED_POINTS\n";
            out_file << "DIMENSIONS " << Constants->dimX_int << " " << Constants->dimY_int << " " << Constants->dimZ_int
                     << "\n";
            out_file << "ASPECT_RATIO " << 1.f / Constants->resolution_factor << " "
                     << 1.f / Constants->resolution_factor
                     << " "
                     << 1.f / Constants->resolution_factor << "\n";
            out_file << "ORIGIN " << Constants->BBminx - Constants->additional_factor << " "
                     << Constants->BBminy - Constants->additional_factor << " "
                     << Constants->BBminz - Constants->additional_factor << "\n";
            out_file << "POINT_DATA " << Constants->dimX_int * Constants->dimY_int * Constants->dimZ_int << "\n";

            out_file << "SCALARS " << "rho " << "double 1\n";
            out_file << "LOOKUP_TABLE default\n";

            for (int k = 0; k < Constants->dimZ_int; k++) {
                for (int j = 0; j < Constants->dimY_int; j++) {
                    for (int i = 0; i < Constants->dimX_int; i++) {
                        out_file << Data->rho(i, j, k) * Constants->Crho << "\n";
                    }
                }
            }

            out_file << "SCALARS " << "rho " << "double 1\n";
            out_file << "LOOKUP_TABLE default\n";

            for (int k = 0; k < Constants->dimZ_int; k++) {
                for (int j = 0; j < Constants->dimY_int; j++) {
                    for (int i = 0; i < Constants->dimX_int; i++) {
                        out_file << Data->p(i, j, k) * Constants->Cpressure << "\n";
                    }
                }
            }


            out_file << "VECTORS " << "U " << "double\n";

            for (int k = 0; k < Constants->dimZ_int; k++) {
                for (int j = 0; j < Constants->dimY_int; j++) {
                    for (int i = 0; i < Constants->dimX_int; i++) {
                        out_file << Data->u(i, j, k).x() * Constants->Cu << " "
                                 << Data->u(i, j, k).y() * Constants->Cu << " "
                                 << Data->u(i, j, k).z() * Constants->Cu << "\n";
                    }
                }
            }

            out_file.close();

        }

    }


};

#endif //INC_24_4_GEO_TO_TNL_OUTPUTER_H
