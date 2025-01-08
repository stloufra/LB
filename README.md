# Lattice Boltzmann Method (LBM) Simulation  

This repository provides a dimploma thesis framework for Lattice Boltzmann Method (LBM) simulation.

## Installation  

### Prerequisites  

Ensure the following software and libraries are installed:  

- **C++ Compiler**: A modern C++17-compatible compiler (e.g., GCC 9+, Clang 10+). For faster simulation preparation C++20 recomended.   
- **CUDA Toolkit**: Required for GPU acceleration (CUDA 11.0 or later recommended). nvcc  
- **CMake**: Version 3.18 or higher.  
- **TNL Library**: Install Template Numerical Library following [TNL Documentation](ttps://tnl-project.gitlab.io/tnl/index.html#installation).
- **CGAL**: Computational Geometry Algorithms Library for preprocess. Tested with version 5.6 .

- **ParaView (Optional)**: For postprocessing and visualization.  

### Build Instructions  

Assuming CUDA capable device availability.

#### 1. Clone the Repository  

```bash
git clone https://github.com/stloufra/LB.git
git checkout thesis
cd LB
```

The framework consists of two main parts: 

1. **2D**
2. **3D**

and 3D is futher divided into:

1. **Preprocessing**: Setting up the simulation domain from .OFF file.  
2. **Simulation**: Running the LBM simulation.

#### 2. Run 2D example

2D is prepared with simple laminar example of flow around cylinder with symetry on one wall.

How to run:
```bash
cd ./2D
make run
```
Results are stored in folder results in .vtk format and can be easily viewed with ParaView.

#### 2. Run 3D example

3D is prepared with LES simulation of backward facing step.

How to run:

1. Preprocess the .OFF file
   
```bash
cd ./3D/preprocess
make run
```
Preprocessed lattice like mesh can be checked in .vtk format file in ./3D/simulation/results.

2. Run the simulation
   
```bash
cd ./3D/simulation
make run
```

Results are stored in folder results in .vtk format and can be easily viewed with ParaView.
