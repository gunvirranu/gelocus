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

static const double PI = 3.141592653589793238462643383279;
static const double RAD_PER_DEG = 0.017453292519943295769;
static const double DEG_PER_RAD = 57.29577951308232087680;
static const double RAD_PER_ARCSEC = PI / (180.0 * 3600);

static inline double rad_to_deg(const double rad) {
    return rad * DEG_PER_RAD;
}

static inline double deg_to_rad(const double deg) {
    return deg * RAD_PER_DEG;
}

#endif  // GELOCUS_COMMON_PRIVATE_H
