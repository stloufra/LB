#include <iostream>
#include <fstream>
#include <string>
#include <list>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include "../traits/LBMTraits.h"
#include "../data/LBMConstants.h"

#pragma once


typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Vector_3 Vector3;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef K::Ray_3 Ray;

typedef CGAL::Polyhedron_3<K> Polyhedron;

typedef CGAL::AABB_polyhedron_triangle_primitive<K, Polyhedron> Primitive;

typedef CGAL::AABB_traits<K, Primitive> Traits;

typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Object_and_primitive_id Object_and_primitive_id;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

typedef CGAL::Bbox_3 Bbox_3;



struct Polytree {
    Polyhedron *poly;
    Tree *tree;
};

// TNL
using DeviceType = typename LBMTraits::DeviceType;
using LBMConstantsPointer = TNL::Pointers::SharedPointer< LBMConstants, DeviceType >;

#ifndef GEOMETRYHANDLEROFF_H
#define GEOMETRYHANDLEROFF_H

#pragma once

extern "C" {

class geometryHandlerOFF {
public:
    geometryHandlerOFF() {
                error = 0;

    }

    void polyhedronFromFile( LBMConstantsPointer Constants, bool verbose) {
        const char *fname = Constants -> file_name.c_str();

        Name = Constants -> file_name;

        Polytree **pptree = &ppt;

        Polyhedron *polyhedron = new Polyhedron;

        std::ifstream in(fname);

        if (verbose) { std::cout << " Reading file " << fname << " " << std::endl; }

        try {
            in >> *polyhedron;
        }
        catch (...) {
            error = 2;
            return;
        }

        Tree *tree = new Tree(polyhedron->facets_begin(), polyhedron->facets_end());

        if (verbose) {
            std::cout << " facets: " << polyhedron->size_of_facets() << std::endl;
            std::cout << " halfedges: " << polyhedron->size_of_halfedges() << std::endl;
            std::cout << " vertices: " << polyhedron->size_of_vertices() << std::endl;
        }


        if (polyhedron->size_of_facets() == 0 ||
            polyhedron->size_of_halfedges() == 0 ||
            polyhedron->size_of_vertices() == 0) {
            error = 1;
            return;
        };

        tree->accelerate_distance_queries();

        ppt = new Polytree;
        (*pptree)->poly = polyhedron;
        (*pptree)->tree = tree;

        error = 0;

    }

    void polyhedronBbox(bool verbose) {
        const Polytree *ptree = ppt;
        d3 *const min = &d3_min_BB;
        d3 *const max = &d3_max_BB;
        Bbox_3 bbox = ptree->tree->bbox();
        *min = {bbox.xmin(), bbox.ymin(), bbox.zmin()};
        *max = {bbox.xmax(), bbox.ymax(), bbox.zmax()};

        if (verbose) {
            std::cout << "BBox of " << Name << " is [" << d3_min_BB.x << ", " << d3_min_BB.y << ", " << d3_min_BB.z
                      << "] to [" << d3_max_BB.x << ", " << d3_max_BB.y << ", " << d3_max_BB.z << "].\n";
        }
    }

    void writeToConstants(LBMConstantsPointer& Constants)
    {
                Constants -> BBminx=d3_min_BB.x;
                Constants -> BBmaxx=d3_max_BB.x;
                Constants -> BBminy=d3_min_BB.y;
                Constants -> BBmaxy=d3_max_BB.y;
                Constants -> BBminz=d3_min_BB.z;
                Constants -> BBmaxz=d3_max_BB.z;

    }

    bool polyhedronInside(d3 pnt_ask, const d3 pnt_ref) {
        const d3 *query = &pnt_ask;
        const d3 *outside_ref = &pnt_ref;

        if (query->x == outside_ref->x && query->y == outside_ref->y && query->z == outside_ref->z) {
            return 0;
        }

        Segment seg(Point(query->x, query->y, query->z),
                    Point(outside_ref->x, outside_ref->y, outside_ref->z));

        const Polytree *ptree = ppt;

        std::list<Object_and_primitive_id> intersections;

        ptree->tree->all_intersections(seg, std::back_inserter(intersections));


        std::vector<Point> points;

        int i = 0;
        for (auto iter = intersections.begin(); iter != intersections.end(); ++iter) {
            i += 1;
            // gets intersection object
            Object_and_primitive_id op = *iter;
            CGAL::Object object = op.first;
            Point point;
            if (CGAL::assign(point, object)) {
                points.push_back(point);
            }
        }

        int n_dist = 0;
        // find how many of the points are distinct
        for (std::vector<Point>::size_type i = 0; i < points.size(); ++i) {
            bool distinct = true;

            for (std::vector<Point>::size_type j = 0; j < i; ++j) {
                Vector3 v = points[i] - points[j];
                distinct = (v.squared_length() > 1e-10);

                if (!distinct) break;
            }
            if (distinct) n_dist += 1;
        }

        return n_dist % 2 == 1;
    }

    ~geometryHandlerOFF() {
        Polytree **pptree = &ppt;
        delete (*pptree)->tree;
        (*pptree)->tree = NULL;
        delete (*pptree)->poly;
        (*pptree)->poly = NULL;
        delete *pptree;
        *pptree = NULL;
    };

    d3 d3_min_BB;
    d3 d3_max_BB;


private:

    int error;
    Polytree *ppt;
    std::string Name;

};

}

#endif