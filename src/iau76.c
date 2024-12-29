/*
 * gelocus -- Earth-centered reference frame transformations library
 *
 * "Scientific code that shouldn't feel like it was written by scientists!"
 *
 * Licensed under the MIT License <http://opensource.org/licenses/MIT>.
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2024 Gunvir Singh Ranu
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "gelocus/iau76.h"

#include <stddef.h>
#include <math.h>

#include "common_private.h"

double lib_gelocus_gmst(const double jc_ut1) {
    static const double cs[] = {
            67310.54841, 876600.0 * 3600 + 8640184.812866, 0.093104, -6.2e-6
    };
    // Horner's method for polynomial evaluation
    const double gmst_secs = cs[0] + jc_ut1 * (cs[1] + jc_ut1 * (cs[2] + jc_ut1 * cs[3]));
    double gmst = fmod(deg_to_rad(gmst_secs / 240), 2 * PI);
    if (gmst < 0) {
        gmst += 2 * PI;
    }
    return gmst;
}

void lib_gelocus_iau76_precession(const double jc, lib_gelocus_Matrix3 * const P)
{
    if (P == NULL)
    {
        return;
    }

    // TODO: explain where #'s came from
    // All in [arcsecond], uses Horner's method
    const double zeta_as = jc * (2306.2181 + jc * (0.30188 + jc * 0.017998));
    const double theta_as = jc * (2004.3109 + jc * (-0.42665 + jc * -0.041833));
    const double z_as = jc * (2306.2181 + jc * (1.09468 + jc * 0.018203));

    // Convert all from [arcsecond] to [rad]
    const double zeta = zeta_as * RAD_PER_ARCSEC;
    const double theta = theta_as * RAD_PER_ARCSEC;
    const double z = z_as * RAD_PER_ARCSEC;

    const double cos_zeta = cos(zeta);
    const double sin_zeta = sin(zeta);
    const double cos_theta = cos(theta);
    const double sin_theta = sin(theta);
    const double cos_z = cos(z);
    const double sin_z = sin(z);

    // TODO: explain rotation choice
    P->row1.x = cos_zeta * cos_theta * cos_z - sin_zeta * sin_z;
    P->row1.y = cos_zeta * cos_theta * sin_z + sin_zeta * cos_z;
    P->row1.z = cos_zeta * sin_theta;
    P->row2.x = -sin_zeta * cos_theta * cos_z - cos_zeta * sin_z;
    P->row2.y = -sin_zeta * cos_theta * sin_z + cos_zeta * cos_z;
    P->row2.z = -sin_zeta * sin_theta;
    P->row3.x = -sin_theta * cos_z;
    P->row3.y = -sin_theta * sin_z;
    P->row3.z = cos_theta;
}

// Moved this ugly boy to the bottom of the file to keep out of the way
const lib_gelocus_IAU80NutationCoeffSet LIB_GELOCUS_IAU80_NUTATION_COEFFS[LIB_GELOCUS_IAU80_NUTATION_TERMS] = {
    {  0,  0,  0,  0,  1, -171996.0, -174.2, 92025.0,  8.9,   1 },
    {  0,  0,  2, -2,  2,  -13187.0,   -1.6,  5736.0, -3.1,   9 },
    {  0,  0,  2,  0,  2,   -2274.0,   -0.2,   977.0, -0.5,  31 },
    {  0,  0,  0,  0,  2,    2062.0,    0.2,  -895.0,  0.5,   2 },
    {  0,  1,  0,  0,  0,    1426.0,   -3.4,    54.0, -0.1,  10 },
    {  1,  0,  0,  0,  0,     712.0,    0.1,    -7.0,  0.0,  32 },
    {  0,  1,  2, -2,  2,    -517.0,    1.2,   224.0, -0.6,  11 },
    {  0,  0,  2,  0,  1,    -386.0,   -0.4,   200.0,  0.0,  33 },
    {  1,  0,  2,  0,  2,    -301.0,    0.0,   129.0, -0.1,  34 },
    {  0, -1,  2, -2,  2,     217.0,   -0.5,   -95.0,  0.3,  12 },
    {  1,  0,  0, -2,  0,    -158.0,    0.0,    -1.0,  0.0,  35 },
    {  0,  0,  2, -2,  1,     129.0,    0.1,   -70.0,  0.0,  13 },
    { -1,  0,  2,  0,  2,     123.0,    0.0,   -53.0,  0.0,  36 },
    {  1,  0,  0,  0,  1,      63.0,    0.1,   -33.0,  0.0,  38 },
    {  0,  0,  0,  2,  0,      63.0,    0.0,    -2.0,  0.0,  37 },
    { -1,  0,  2,  2,  2,     -59.0,    0.0,    26.0,  0.0,  40 },
    { -1,  0,  0,  0,  1,     -58.0,   -0.1,    32.0,  0.0,  39 },
    {  1,  0,  2,  0,  1,     -51.0,    0.0,    27.0,  0.0,  41 },
    {  2,  0,  0, -2,  0,      48.0,    0.0,     1.0,  0.0,  14 },
    { -2,  0,  2,  0,  1,      46.0,    0.0,   -24.0,  0.0,   3 },
    {  0,  0,  2,  2,  2,     -38.0,    0.0,    16.0,  0.0,  42 },
    {  2,  0,  2,  0,  2,     -31.0,    0.0,    13.0,  0.0,  45 },
    {  2,  0,  0,  0,  0,      29.0,    0.0,    -1.0,  0.0,  43 },
    {  1,  0,  2, -2,  2,      29.0,    0.0,   -12.0,  0.0,  44 },
    {  0,  0,  2,  0,  0,      26.0,    0.0,    -1.0,  0.0,  46 },
    {  0,  0,  2, -2,  0,     -22.0,    0.0,     0.0,  0.0,  15 },
    { -1,  0,  2,  0,  1,      21.0,    0.0,   -10.0,  0.0,  47 },
    {  0,  2,  0,  0,  0,      17.0,   -0.1,     0.0,  0.0,  16 },
    {  0,  2,  2, -2,  2,     -16.0,    0.1,     7.0,  0.0,  18 },
    { -1,  0,  0,  2,  1,      16.0,    0.0,    -8.0,  0.0,  48 },
    {  0,  1,  0,  0,  1,     -15.0,    0.0,     9.0,  0.0,  17 },
    {  1,  0,  0, -2,  1,     -13.0,    0.0,     7.0,  0.0,  49 },
    {  0, -1,  0,  0,  1,     -12.0,    0.0,     6.0,  0.0,  19 },
    {  2,  0, -2,  0,  0,      11.0,    0.0,     0.0,  0.0,   4 },
    { -1,  0,  2,  2,  1,     -10.0,    0.0,     5.0,  0.0,  50 },
    {  1,  0,  2,  2,  2,      -8.0,    0.0,     3.0,  0.0,  54 },
    {  0, -1,  2,  0,  2,      -7.0,    0.0,     3.0,  0.0,  53 },
    {  0,  0,  2,  2,  1,      -7.0,    0.0,     3.0,  0.0,  58 },
    {  1,  1,  0, -2,  0,      -7.0,    0.0,     0.0,  0.0,  51 },
    {  0,  1,  2,  0,  2,       7.0,    0.0,    -3.0,  0.0,  52 },
    { -2,  0,  0,  2,  1,      -6.0,    0.0,     3.0,  0.0,  20 },
    {  0,  0,  0,  2,  1,      -6.0,    0.0,     3.0,  0.0,  57 },
    {  2,  0,  2, -2,  2,       6.0,    0.0,    -3.0,  0.0,  56 },
    {  1,  0,  0,  2,  0,       6.0,    0.0,     0.0,  0.0,  55 },
    {  1,  0,  2, -2,  1,       6.0,    0.0,    -3.0,  0.0,  58 },
    {  0,  0,  0, -2,  1,      -5.0,    0.0,     3.0,  0.0,  60 },
    {  0, -1,  2, -2,  1,      -5.0,    0.0,     3.0,  0.0,  21 },
    {  2,  0,  2,  0,  1,      -5.0,    0.0,     3.0,  0.0,  62 },
    {  1, -1,  0,  0,  0,       5.0,    0.0,     0.0,  0.0,  61 },
    {  1,  0,  0, -1,  0,      -4.0,    0.0,     0.0,  0.0,  24 },
    {  0,  0,  0,  1,  0,      -4.0,    0.0,     0.0,  0.0,  65 },
    {  0,  1,  0, -2,  0,      -4.0,    0.0,     0.0,  0.0,  63 },
    {  1,  0, -2,  0,  0,       4.0,    0.0,     0.0,  0.0,  64 },
    {  2,  0,  0, -2,  1,       4.0,    0.0,    -2.0,  0.0,  22 },
    {  0,  1,  2, -2,  1,       4.0,    0.0,    -2.0,  0.0,  23 },
    {  1,  1,  0,  0,  0,      -3.0,    0.0,     0.0,  0.0,  66 },
    {  1, -1,  0, -1,  0,      -3.0,    0.0,     0.0,  0.0,   6 },
    { -1, -1,  2,  2,  2,      -3.0,    0.0,     1.0,  0.0,  69 },
    {  0, -1,  2,  2,  2,      -3.0,    0.0,     1.0,  0.0,  72 },
    {  1, -1,  2,  0,  2,      -3.0,    0.0,     1.0,  0.0,  68 },
    {  3,  0,  2,  0,  2,      -3.0,    0.0,     1.0,  0.0,  71 },
    { -2,  0,  2,  0,  2,      -3.0,    0.0,     1.0,  0.0,   5 },
    {  1,  0,  2,  0,  0,       3.0,    0.0,     0.0,  0.0,  67 },
    { -1,  0,  2,  4,  2,      -2.0,    0.0,     1.0,  0.0,  82 },
    {  1,  0,  0,  0,  2,      -2.0,    0.0,     1.0,  0.0,  76 },
    { -1,  0,  2, -2,  1,      -2.0,    0.0,     1.0,  0.0,  74 },
    {  0, -2,  2, -2,  1,      -2.0,    0.0,     1.0,  0.0,   7 },
    { -2,  0,  0,  0,  1,      -2.0,    0.0,     1.0,  0.0,  70 },
    {  2,  0,  0,  0,  1,       2.0,    0.0,    -1.0,  0.0,  75 },
    {  3,  0,  0,  0,  0,       2.0,    0.0,     0.0,  0.0,  77 },
    {  1,  1,  2,  0,  2,       2.0,    0.0,    -1.0,  0.0,  73 },
    {  0,  0,  2,  1,  2,       2.0,    0.0,    -1.0,  0.0,  78 },
    {  1,  0,  0,  2,  1,      -1.0,    0.0,     0.0,  0.0,  91 },
    {  1,  0,  2,  2,  1,      -1.0,    0.0,     1.0,  0.0,  85 },
    {  1,  1,  0, -2,  1,      -1.0,    0.0,     0.0,  0.0, 102 },
    {  0,  1,  0,  2,  0,      -1.0,    0.0,     0.0,  0.0,  99 },
    {  0,  1,  2, -2,  0,      -1.0,    0.0,     0.0,  0.0,  30 },
    {  0,  1, -2,  2,  0,      -1.0,    0.0,     0.0,  0.0,  27 },
    {  1,  0, -2,  2,  0,      -1.0,    0.0,     0.0,  0.0, 103 },
    {  1,  0, -2, -2,  0,      -1.0,    0.0,     0.0,  0.0, 100 },
    {  1,  0,  2, -2,  0,      -1.0,    0.0,     0.0,  0.0,  94 },
    {  1,  0,  0, -4,  0,      -1.0,    0.0,     0.0,  0.0,  80 },
    {  2,  0,  0, -4,  0,      -1.0,    0.0,     0.0,  0.0,  83 },
    {  0,  0,  2,  4,  2,      -1.0,    0.0,     0.0,  0.0, 105 },
    {  0,  0,  2, -1,  2,      -1.0,    0.0,     0.0,  0.0,  98 },
    { -2,  0,  2,  4,  2,      -1.0,    0.0,     1.0,  0.0,  86 },
    {  2,  0,  2,  2,  2,      -1.0,    0.0,     0.0,  0.0,  90 },
    {  0, -1,  2,  0,  1,      -1.0,    0.0,     0.0,  0.0, 101 },
    {  0,  0, -2,  0,  1,      -1.0,    0.0,     0.0,  0.0,  97 },
    {  0,  0,  4, -2,  2,       1.0,    0.0,     0.0,  0.0,  92 },
    {  0,  1,  0,  0,  2,       1.0,    0.0,     0.0,  0.0,  28 },
    {  1,  1,  2, -2,  2,       1.0,    0.0,    -1.0,  0.0,  84 },
    {  3,  0,  2, -2,  2,       1.0,    0.0,     0.0,  0.0,  93 },
    { -2,  0,  2,  2,  2,       1.0,    0.0,    -1.0,  0.0,  81 },
    { -1,  0,  0,  0,  2,       1.0,    0.0,    -1.0,  0.0,  79 },
    {  0,  0, -2,  2,  1,       1.0,    0.0,     0.0,  0.0,  26 },
    {  0,  1,  2,  0,  1,       1.0,    0.0,     0.0,  0.0,  95 },
    { -1,  0,  4,  0,  2,       1.0,    0.0,     0.0,  0.0,  87 },
    {  2,  1,  0, -2,  0,       1.0,    0.0,     0.0,  0.0,  25 },
    {  2,  0,  0,  2,  0,       1.0,    0.0,     0.0,  0.0, 104 },
    {  2,  0,  2, -2,  1,       1.0,    0.0,    -1.0,  0.0,  89 },
    {  2,  0, -2,  0,  1,       1.0,    0.0,     0.0,  0.0,   8 },
    {  1, -1,  0, -2,  0,       1.0,    0.0,     0.0,  0.0,  88 },
    { -1,  0,  0,  1,  1,       1.0,    0.0,     0.0,  0.0,  29 },
    { -1, -1,  0,  2,  1,       1.0,    0.0,     0.0,  0.0,  96 },
    {  0,  1,  0,  1,  0,       1.0,    0.0,     0.0,  0.0, 106 },
};
