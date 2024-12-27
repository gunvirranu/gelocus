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

#include "gelocus/gelocus.h"

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
