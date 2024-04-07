#ifndef OBJ_RECTANGLE_H
#define OBJ_RECTANGLE_H
#include <cmath> 
#include "Obj_template.h"

#pragma once

class Obj_rectangle : public Obj_template
{
public:
    Obj_rectangle()
    {
        x_min = 0;
        x_max = 0;
        y_min = 0;
        y_max = 0;
    
    }

    Obj_rectangle(double x_min_in , double x_max_in , double y_min_in, double y_max_in)
    {
        x_min = x_min_in;
        x_max = x_max_in;
        y_min = y_min_in;
        y_max = y_max_in;
        
    }

    virtual bool is_inside(double x_b, double y_b) const
    {
        bool vnitr = 0;   
        if(x_b >= x_min && x_b <= x_max && y_b >= y_min && y_b <= y_max)
            vnitr=1;
        else 
            vnitr=0; 

        return vnitr;

    }


    virtual ~Obj_rectangle() {}
    

private:
    double x_min;
    double x_max;
    double y_min;
    double y_max;
};

#endif