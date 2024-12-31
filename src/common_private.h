/**
 * gelocus -- Earth-centered reference frame transformations library
 * https://github.com/gunvirranu/gelocus
 *
 * SPDX-License-Identifier: MIT
 * Copyright (c) 2024 Gunvir Singh Ranu
 */

#ifndef GELOCUS_COMMON_PRIVATE_H
#define GELOCUS_COMMON_PRIVATE_H

#define UNUSED(x)   ((void) x)

static const double PI = 3.14159265358979323846264338328;   ///< [-] Our favourite mathy boi
static const double RAD_TO_DEG = 57.29577951308232087680;   ///< [deg/rad] 180 / PI
static const double DEG_TO_RAD = 0.017453292519943295769;   ///< [rad/deg] PI / 180

static const double DEG_TO_ARCMIN    = 60.0;                ///< [arcmin/deg] Definition of minute of arc
static const double DEG_TO_ARCSEC    = 3600.0;              ///< [arcsec/deg] DEG_TO_ARCMIN * ARCMIN_TO_ARCSEC
static const double ARCMIN_TO_ARCSEC = 60.0;                ///< [arcsec/arcmin] Definition of second of arc

static const double ARCMIN_TO_RAD = 2.908882086657215e-4;   ///< [rad/arcmin] DEG_TO_RAD / DEG_TO_ARCMIN
static const double ARCSEC_TO_RAD = 4.848136811095359e-6;   ///< [rad/arcsec] DEG_TO_RAD / DEG_TO_ARCSEC
static const double MAS_TO_RAD    = 4.848136811095359e-9;   ///< [rad/mas] 1e-3 * ARCSEC_TO_RAD
static const double P1MAS_TO_RAD  = 4.848136811095359e-10;  ///< [rad/(0.1 mas)] 1e-4 * ARCSEC_TO_RAD

static inline double rad_to_deg(const double rad) {
    return rad * RAD_TO_DEG;
}

static inline double deg_to_rad(const double deg) {
    return deg * DEG_TO_RAD;
}

#endif  // GELOCUS_COMMON_PRIVATE_H
