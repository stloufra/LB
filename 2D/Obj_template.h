#ifndef OBJ_TEMPLATE_H
#define OBJ_TEMPLATE_H

#pragma once

class Obj_template
{
public:
    virtual bool is_inside(double x_b,double y_b) const = 0;
    
};

#endif