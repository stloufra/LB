//
// Created by stloufra on 10/30/23.
//

#ifndef D3Q15_H
#define D3Q15_H

#include "../../../traits/LBMTraits.h"

struct D3Q15
{
    D3Q15() = default;

    using RealType = LBMTraits::RealType;

    static constexpr const int  numberOfDiscreteVelocities = 15;
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

        // weight is 2/9 for 0
        ooo = 0,

        // weight is 2/27 for (1-6)
        eoo = 1,
        woo = 2,
        oob = 3,
        oof = 4,
        oso = 5,
        ono = 6,


        // weight is 1/72 for (7-14)
        wnb = 7,
        esf = 8,
        wsf = 9,
        enb = 10,
        esb = 11,
        wnf = 12,
        wsb = 13,
        enf = 14

    };

    const int c[15][3] = {
            {0,  0,  0},

            {1,  0,  0},
            {-1, 0,  0},
            {0,  0,  -1},
            {0,  0,  1},
            {0,  -1, 0},
            {0,  1,  0},

            {-1, 1,  -1},
            {1,  -1, 1},
            {-1, -1, 1},
            {1,  1,  -1},
            {1,  -1, -1},
            {-1, 1,  1},
            {-1, -1, -1},
            {1,  1,  1}};

    const RealType weight[15] = {   2.f / 9.f, // 0

                                    1.f / 9.f, 1.f / 9.f, 1.f / 9.f, 1.f / 9.f, 1.f / 9.f, 1.f / 9.f, // 1-6

                                    1.f / 72.f, 1.f / 72.f, 1.f / 72.f, 1.f / 72.f, 1.f / 72.f,
                                    1.f / 72.f, 1.f / 72.f, 1.f / 72.f}; // 7 - 14

    const int c_rev[15] = {   0, 2, 1, 4, 3, 6,
                              5, 8, 7, 10, 9, 12, 11,
                              14, 13}; // (reverse velocities are i+1 for i odd and i-1 for i even)


};

#endif //SIMULATION_D3Q15_H
