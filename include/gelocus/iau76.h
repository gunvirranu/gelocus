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

#ifndef GELOCUS_IAU76_H
#define GELOCUS_IAU76_H

#include <stddef.h>

#include "gelocus.h"

// IAU-76 / FK5 Theory

#define LIB_GELOCUS_IAU80_NUTATION_TERMS  106U

/// Coefficient set for the IAU-1980 nutation theory
///
/// @todo: detail exactly which model and from where
typedef struct {
    int l;          ///< [-] l  coefficient
    int l1;         ///< [-] l1 coefficient
    int f;          ///< [-] f  coefficient
    int d;          ///< [-] d  coefficient
    int omg;        ///< [-] omega coefficient
    double sp;      ///< [0.1 mas] Longitude sine coefficient
    double spt;     ///< [0.1 mas / julian century] Transient longitude sine coefficient
    double ce;      ///< [0.1 mas] Obliquity cosine coefficient
    double cet;     ///< [0.1 mas / julian century] Transient obliquity cosine coefficients
    size_t idx;     ///< Original index of coefficient set
} lib_gelocus_IAU80NutationCoeffSet;

/// IAU-1980 Nutation theory coefficients
///
/// Ordered by largest coefficients first for summation. Original index is last term `idx`.
extern const lib_gelocus_IAU80NutationCoeffSet LIB_GELOCUS_IAU80_NUTATION_COEFFS[LIB_GELOCUS_IAU80_NUTATION_TERMS];

/// Greenwich Mean Sidereal Time in UT1
double lib_gelocus_gmst(double jc_ut1);

/// IAU-1976 Precession theory for MOD to J2000
void lib_gelocus_iau76_precession(double jc_ut1, lib_gelocus_Matrix3 * P);

#endif  // GELOCUS_IAU76_H
