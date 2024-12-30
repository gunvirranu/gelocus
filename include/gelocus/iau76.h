/**
 * gelocus -- Earth-centered reference frame transformations library
 * https://github.com/gunvirranu/gelocus
 *
 * SPDX-License-Identifier: MIT
 * Copyright (c) 2024 Gunvir Singh Ranu
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
