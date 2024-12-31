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
#include <stdint.h>

#include "gelocus.h"

// IAU-76 / FK5 Theory

#define LIB_GELOCUS_IAU80_NUTATION_TERMS  106

/// Coefficient set for the IAU-1980 nutation theory
///
/// @todo: detail exactly which model and from where
typedef struct {
    int8_t l;       ///< [-] l  coefficient
    int8_t l1;      ///< [-] l1 coefficient
    int8_t f;       ///< [-] f  coefficient
    int8_t d;       ///< [-] d  coefficient
    int8_t omg;     ///< [-] omega coefficient
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

/// FK5 Polar Motion theory for ECEF to PEF
void lib_gelocus_fk5_polar_motion(
    lib_gelocus_EOPData eop,    ///< Earth Orientation Parameters
    lib_gelocus_Matrix3 * PM    ///< Polar motion rotation matrix (ignored if NULL)
);

/// SGP4 TEME to PEF transformation
void lib_gelocus_teme_to_pef(
    double jc_ut1,              ///< [julian century] Timestamp in UT1
    lib_gelocus_Matrix3 * S     ///< TEME to PEF rotation matrix (ignored if NULL)
);

/// IAU-1976 Nutation theory for PEF to TOD
void lib_gelocus_iau76_sidereal(
    double jc_ut1,              ///< [julian century] Timestamp in UT1
    double mean_eps,            ///< @todo
    double omega,               ///< @todo
    double delta_psi,           ///< @todo
    lib_gelocus_Matrix3 * S     ///< Sidereal rotation matrix (ignored if NULL)
);

/// IAU-1980 Nutation theory for TOD to MOD
///
/// @note The pointer outputs are optional, as in they are ignored if NULL.
///       You are welcome to only consume those which you desire.
void lib_gelocus_iau80_nutation(
    double jc_ut1,              ///< [julian century] Timestamp in UT1
    lib_gelocus_EOPData eop,    ///< Earth Orientation Parameters
    lib_gelocus_Matrix3 * N,    ///< Nutation rotation matrix (ignored if NULL)
    // TODO: consider consolidating into a struct
    double * mean_eps,          ///< @todo (ignored if NULL)
    double * omega,             ///< @todo (ignored if NULL)
    double * delta_psi          ///< @todo (ignored if NULL)
);

/// IAU-1976 Precession theory for MOD to J2000
void lib_gelocus_iau76_precession(double jc_ut1, lib_gelocus_Matrix3 * P);

#endif  // GELOCUS_IAU76_H
