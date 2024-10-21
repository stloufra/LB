//
// Created by stloufra on 10/30/23.
//

#ifndef D3Q19_H
#define D3Q19_H

#include "../../../traits/LBMTraits.h"

struct D3Q19
{
    D3Q19() = default;

    using RealType = LBMTraits::RealType;

    static constexpr const int  numberOfDiscreteVelocities = 19;
    enum {
        //velocities are denoted by 4 cardinal directions:
        // "n" - north, "s" - south, "e" - east and "w" - west
        // in addition with "f" - forward, "b" - backward direction,
        // "o" - non-part of vector.
        //
        //           +N2
        //            ↑   -B3
        //            | ⋰
        //  -W1 <-----O-----> +E1
        //          ⋰ |
        //      +F3   ↓
        //           -S2
        //
        //            +y
        //            ↑    -z
        //            | ⋰
        //   -x <-----O-----> +x
        //          ⋰ |
        //      +z    ↓
        //            -y
        //
        // index precedence is:
        // 1: {w, o, e}
        // 2: {s, o, n}
        // 3: {b, o, f}

        // weight is 8/27 for 0
        ooo = 0,

        // weight is 2/27 for (1-6)
        eoo = 1,
        woo = 2,
        oob = 3,
        oof = 4,
        oso = 5,
        ono = 6,

        // weight is 1/54 for (7-18)
        eob = 7,
        wof = 8,
        eof = 9,
        wob = 10,
        wso = 11,
        eno = 12,
        onb = 13,
        osf = 14,
        wno = 15,
        eso = 16,
        onf = 17,
        osb = 18,



    };

    const int c[19][3] = {
            {0,  0,  0},

            {1,  0,  0},
            {-1, 0,  0},
            {0,  0,  -1},
            {0,  0,  1},
            {0,  -1, 0},
            {0,  1,  0},

            {1,  0,  -1},
            {-1, 0,  1},
            {1,  0,  1},
            {-1, 0,  -1},
            {-1, -1, 0},
            {1,  1,  0},
            {0,  1,  -1},
            {0,  -1, 1},
            {-1, 1,  0},
            {1,  -1, 0},
            {0,  1,  1},
            {0,  -1, -1}};



    const RealType weight[19] = {1.f / 3.f, // 0
                1.f / 18.f, 1.f / 18.f, 1.f / 18.f, 1.f / 18.f, 1.f / 18.f, 1.f / 18.f, // 1-6
                1.f / 36.f, 1.f / 36.f, 1.f / 36.f, 1.f / 36.f, 1.f / 36.f, 1.f / 36.f,
                1.f / 36.f, 1.f / 36.f, 1.f / 36.f, 1.f / 36.f, 1.f / 36.f, 1.f / 36.f};// 7- 18


    const int c_rev[19] = {0, 2, 1, 4, 3, 6,
                                      5, 8, 7, 10, 9, 12, 11,
                                      14, 13, 16, 15, 18, 17}; // (reverse velocities are i+1 for i odd and i-1 for i even)

    const int c_sym_x[27] = {
        0,   // 0 -> {0, 0, 0} stays the same

        2,   // 1 -> {1, 0, 0} becomes {-1, 0, 0} (index 2)
        1,   // 2 -> {-1, 0, 0} becomes {1, 0, 0} (index 1)
        3,   // 3 -> {0, 0, -1} stays the same
        4,   // 4 -> {0, 0, 1} stays the same
        5,   // 5 -> {0, -1, 0} stays the same
        6,   // 6 -> {0, 1, 0} stays the same

        10,  // 7 -> {1, 0, -1} becomes {-1, 0, -1} (index 10)
        9,   // 8 -> {-1, 0, 1} becomes {1, 0, 1} (index 9)
        8,   // 9 -> {1, 0, 1} becomes {-1, 0, 1} (index 8)
        7,   // 10 -> {-1, 0, -1} becomes {1, 0, -1} (index 7)
        11,  // 11 -> {-1, -1, 0} becomes {1, -1, 0} (index 16)
        12,  // 12 -> {1, 1, 0} becomes {-1, 1, 0} (index 15)
        13,  // 13 -> {0, 1, -1} stays the same
        14,  // 14 -> {0, -1, 1} stays the same
        12,  // 15 -> {-1, 1, 0} becomes {1, 1, 0} (index 12)
        11,  // 16 -> {1, -1, 0} becomes {-1, -1, 0} (index 11)
        17,  // 17 -> {0, 1, 1} stays the same
        18,  // 18 -> {0, -1, -1} stays the same
    };

    const int c_sym_y[27] = {
        0,   // 0 -> {0, 0, 0} stays the same

        1,   // 1 -> {1, 0, 0} stays the same
        2,   // 2 -> {-1, 0, 0} stays the same
        3,   // 3 -> {0, 0, -1} stays the same
        4,   // 4 -> {0, 0, 1} stays the same
        6,   // 5 -> {0, -1, 0} becomes {0, 1, 0} (index 6)
        5,   // 6 -> {0, 1, 0} becomes {0, -1, 0} (index 5)

        7,   // 7 -> {1, 0, -1} stays the same
        8,   // 8 -> {-1, 0, 1} stays the same
        9,   // 9 -> {1, 0, 1} stays the same
        10,  // 10 -> {-1, 0, -1} stays the same
        16,  // 11 -> {-1, -1, 0} becomes {-1, 1, 0} (index 15)
        15,  // 12 -> {1, 1, 0} becomes {1, -1, 0} (index 16)
        18,  // 13 -> {0, 1, -1} becomes {0, -1, -1} (index 18)
        17,  // 14 -> {0, -1, 1} becomes {0, 1, 1} (index 17)
        11,  // 15 -> {-1, 1, 0} becomes {-1, -1, 0} (index 11)
        12,  // 16 -> {1, -1, 0} becomes {1, 1, 0} (index 12)
        14,  // 17 -> {0, 1, 1} becomes {0, -1, 1} (index 14)
        13,  // 18 -> {0, -1, -1} becomes {0, 1, -1} (index 13)
    };

    const int c_sym_z[27] = {
        0,   // 0 -> {0, 0, 0} stays the same

        1,   // 1 -> {1, 0, 0} stays the same
        2,   // 2 -> {-1, 0, 0} stays the same
        4,   // 3 -> {0, 0, -1} becomes {0, 0, 1} (index 4)
        3,   // 4 -> {0, 0, 1} becomes {0, 0, -1} (index 3)
        5,   // 5 -> {0, -1, 0} stays the same
        6,   // 6 -> {0, 1, 0} stays the same

        9,   // 7 -> {1, 0, -1} becomes {1, 0, 1} (index 9)
        10,  // 8 -> {-1, 0, 1} becomes {-1, 0, -1} (index 10)
        7,   // 9 -> {1, 0, 1} becomes {1, 0, -1} (index 7)
        8,   // 10 -> {-1, 0, -1} becomes {-1, 0, 1} (index 8)
        11,  // 11 -> {-1, -1, 0} stays the same
        12,  // 12 -> {1, 1, 0} stays the same
        18,  // 13 -> {0, 1, -1} becomes {0, 1, 1} (index 17)
        17,  // 14 -> {0, -1, 1} becomes {0, -1, -1} (index 18)
        15,  // 15 -> {-1, 1, 0} stays the same
        16,  // 16 -> {1, -1, 0} stays the same
        14,  // 17 -> {0, 1, 1} becomes {0, 1, -1} (index 13)
        13,  // 18 -> {0, -1, -1} becomes {0, -1, 1} (index 14)
    };

};

#endif //SIMULATION_D3Q27DATA_H
