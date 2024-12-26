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

const double LIB_GELOCUS_EPOCH_J2000_JD = 2451545.0;
const double LIB_GELOCUS_EPOCH_J1900_JD = 2415021.0;
const double LIB_GELOCUS_EPOCH_GPS_JD   = 2444244.5;

const double LIB_GELOCUS_EPOCH_J2000_JC = 0.0;

const double LIB_GELOCUS_DELTA_JD_PER_DAY = 1.0;
const double LIB_GELOCUS_DELTA_JC_PER_DAY = 1.0 / 36525;
const double LIB_GELOCUS_DELTA_DAY_PER_JC = 36525;

const lib_gelocus_Matrix3 LIB_GELOCUS_MATRIX_IDENTITY = {
    .row1 = { 1, 0, 0 },
    .row2 = { 0, 1, 0 },
    .row3 = { 0, 0, 1 }
};

double lib_gelocus_vec_norm(const lib_gelocus_Vec3 x)
{
    return sqrt(lib_gelocus_dot_product(x, x));
}

double lib_gelocus_dot_product(const lib_gelocus_Vec3 a, const lib_gelocus_Vec3 b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

lib_gelocus_Vec3 lib_gelocus_multiply_matrix(
    const lib_gelocus_Matrix3 mat,
    const lib_gelocus_Vec3 vec
) {
    lib_gelocus_Vec3 out = { 0 };
    out.x = lib_gelocus_dot_product(mat.row1, vec);
    out.y = lib_gelocus_dot_product(mat.row2, vec);
    out.z = lib_gelocus_dot_product(mat.row3, vec);
    return out;
}

double lib_gelocus_jd_to_jc(const double jd)
{
    return (jd - LIB_GELOCUS_EPOCH_J2000_JD) / LIB_GELOCUS_DELTA_DAY_PER_JC;
}

double lib_gelocus_jd_frac_to_jc(const double jd, const double jd_frac)
{
    const double delta_jd = (jd - LIB_GELOCUS_EPOCH_J2000_JD) + jd_frac;
    return delta_jd / LIB_GELOCUS_DELTA_DAY_PER_JC;
}

bool lib_gelocus_apply_transformation(
    const lib_gelocus_Transformation trans,
    const lib_gelocus_StateVector sv,
    lib_gelocus_StateVector * const out
) {
    const bool is_bad_ptr = (out == NULL);
    const bool is_bad_frame = (trans.from != sv.frame);
    const bool is_bad_time = (trans.jc != sv.jc);

    if (is_bad_ptr || is_bad_frame || is_bad_time)
    {
        return false;
    }

    // Apply linear transformation
    out->pos = lib_gelocus_multiply_matrix(trans.matrix, sv.pos);
    /// FIXME: also do out->vel

    // Copy over stuff that doesn't change
    out->jc = sv.jc;
    out->frame = trans.to;
    out->eop = trans.eop;
    return true;
}
