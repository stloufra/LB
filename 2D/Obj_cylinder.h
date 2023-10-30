#ifndef OBJ_CYLINDER_H
#define OBJ_CYLINDER_H
#include <cmath>
#include "Obj_template.h"

#pragma once


class Obj_cylinder : public Obj_template
{
    public:
        Obj_cylinder()
        {
            r = 0.0;
            x = 0.0;
            y = 0.0;
        }

        Obj_cylinder(double r_in, double x_in, double y_in )
        {
            r = r_in;
            x = x_in;
            y = y_in;
            
        }

        virtual bool is_inside(double x_b, double y_b) const
        {
            bool vnitr = 0;
            if(std::sqrt((x_b-x)*(x_b-x) + (y_b-y)*(y_b-y)) < r)
                vnitr=1;
            else 
                vnitr=0; 

            return vnitr;
        }

        virtual ~Obj_cylinder() {}

    private:
        double r;
        double x;
        double y;
        double const pi = 3.1415;
};

#endif