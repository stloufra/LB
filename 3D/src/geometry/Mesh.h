#ifndef MESH_H
#define MESH_H

#include <vector>
#include <cassert>

//typedef struct {double x,y,z;} d3; // declare in a geometryHandler.h, move to all definitions after
#pragma once


/*
       Mesh stores the orthogonal, evenly-spaced grid of points.
           Point values represent:
           - 0 - solid(wall)
           - 1 - fluid
           - 2 - inlet
           - 3 - outlet
           - 4 - moving wall

           Dimensions:
           - x - columns
           - y - rows
           - z - slices

       For initialization, the number of points in each dimension is given (x, y, z).
       Additional_factor is the number of additional points to be added to the grid in each dimension to each side of the x, y and z dimensions given (minimum one).
       By default, operators do not add additional points to the grid access.
       Scale factor is just saved for visualization .
*/
class Mesh {

public:

    Mesh() = delete;


    Mesh(int x, int y, int z, int additional_factor, int scale_factor) {

        add_f = additional_factor;
        scale = static_cast<double> (scale_factor);

        cols_pr = static_cast<int>((x + add_f * 2));
        rows_pr = static_cast<int>((y + add_f + 2));
        slices_pr = static_cast<int>((z + add_f * 2));

        data.resize(rows_pr * cols_pr * slices_pr);
        cols = x;
        rows = y;
        slices = z;
    }


    int &operator()(int i, int j, int k) {
        assert(j >= 0 && j < rows && i >= 0 && i < cols && k >= 0 && k < slices);

        //assert(j >= -additional_factor && j < rows_pr-additional_factor && i >= -additional_factor && i < cols_pr-additional_factor &&  k >= -additional_factor && k < slices_pr-additional_factor );

        return data[(cols_pr * rows_pr) * (add_f + k) + (add_f + j) * cols_pr + (i + add_f)];
    }

    const int &operator()(int i, int j, int k) const {
        assert(j >= 0 && j < rows && i >= 0 && i < cols && k >= 0 && k < slices);

        return data[(cols_pr * rows_pr) * (add_f + k) + (add_f + j) * cols_pr + (i + add_f)];
    }


    ~Mesh() {}

    int rows;
    int cols;
    int slices;
    int add_f;
    double scale;


private:
    std::vector<int> data;
    int rows_pr;
    int cols_pr;
    int slices_pr;


};

#endif