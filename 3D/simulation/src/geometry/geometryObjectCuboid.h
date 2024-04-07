
#ifndef GEOMETRYOBJECTCUBOID_H
#define GEOMETRYOBJECTCUBOID_H
#include <omp.h>
#include "../traits/LBMTraits.h"



class geometryObjectCuboid {

    public:
    geometryObjectCuboid(d3 X1, d3 X2, d3 X3, int id){

        this->X1 = X1;
        this->X2 = X2;
        this->X3 = X3;
        
        this->id = id;
    }
    
    bool isInside(d3 X_ask){

        x_min = std::min(X1.x, std::min(X2.x, X3.x));
        x_max = std::max(X1.x, std::max(X2.x, X3.x));

        y_min = std::min(X1.y, std::min(X2.y, X3.y));
        y_max = std::max(X1.y, std::max(X2.y, X3.y));

        z_min = std::min(X1.z, std::min(X2.z, X3.z));
        z_max = std::max(X1.z, std::max(X2.z, X3.z));

        if(X_ask.x >= x_min && X_ask.x <= x_max && X_ask.y >= y_min && X_ask.y <= y_max && X_ask.z >= z_min && X_ask.z <= z_max){
            return true;
        }
        else{
            return false;
        }
    }

    int id; //should be unique negative integer

    private:

    d3 X1;
    d3 X2;
    d3 X3;

    double x_min, y_min, z_min;
    double x_max, y_max, z_max;


};


#endif //GEOMETRYOBJECTCUBOID_H
