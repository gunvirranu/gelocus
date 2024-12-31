/**
 * gelocus -- Earth-centered reference frame transformations library
 * https://github.com/gunvirranu/gelocus
 *
 * SPDX-License-Identifier: MIT
 * Copyright (c) 2024 Gunvir Singh Ranu
 */

#include "gelocus/iau76.h"

#include <stddef.h>
#include <math.h>

#include "common_private.h"

double lib_gelocus_gmst(const double jc_ut1)
{
    static const double cs[] = {
            67310.54841, 876600.0 * 3600 + 8640184.812866, 0.093104, -6.2e-6
    };

    // Horner's method for polynomial evaluation
    const double gmst_secs = cs[0] + jc_ut1 * (cs[1] + jc_ut1 * (cs[2] + jc_ut1 * cs[3]));

    double gmst = fmod(deg_to_rad(gmst_secs / 240), 2 * PI);
    if (gmst < 0)
    {
        gmst += 2 * PI;
    }
    return gmst;
}

void lib_gelocus_iau76_sidereal(
    const double jc_ut1,
    const double mean_eps,
    const double omega,
    const double delta_psi,
    lib_gelocus_Matrix3 * const S
) {
    if (S == NULL)
    {
        return;
    }

    const double gmst = lib_gelocus_gmst(jc_ut1);
    double ast = gmst + (delta_psi * cos(mean_eps));

    // TODO: explain this pls pls, also test JD -> JC maybe
    static const double JD_AST_DELTA = 2450449.5;
    if (jc_ut1 > lib_gelocus_jd_to_jc(JD_AST_DELTA))
    {
        ast += ARCSEC_TO_RAD * 0.002640 * sin(omega);
        ast += ARCSEC_TO_RAD * 0.000063 * sin(2 * omega);
    }
    ast = fmod(ast, 2 * PI);

    const double sin_ast = sin(ast);
    const double cos_ast = cos(ast);

    // TODO: explain rotation matrix choice
    S->row1.x = cos_ast;  S->row1.y = -sin_ast;  S->row1.z = 0.0;
    S->row2.x = sin_ast;  S->row2.y =  cos_ast;  S->row2.z = 0.0;
    S->row3.x = 0.0    ;  S->row3.y =  0.0    ;  S->row3.z = 1.0;
}

void lib_gelocus_iau76_precession(const double jc, lib_gelocus_Matrix3 * const P)
{
    if (P == NULL)
    {
        return;
    }

    // TODO: explain where #'s came from
    // All in [arcsecond], uses Horner's method
    const double zeta_arcsec = jc * (2306.2181 + jc * (0.30188 + jc * 0.017998));
    const double theta_arcsec = jc * (2004.3109 + jc * (-0.42665 + jc * -0.041833));
    const double z_arcsec = jc * (2306.2181 + jc * (1.09468 + jc * 0.018203));

    // Convert all from [arcsecond] to [rad]
    const double zeta = zeta_arcsec * ARCSEC_TO_RAD;
    const double theta = theta_arcsec * ARCSEC_TO_RAD;
    const double z = z_arcsec * ARCSEC_TO_RAD;

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

void lib_gelocus_iau80_nutation(
    const double jc,
    const lib_gelocus_EOPData eop,
    lib_gelocus_Matrix3 * const N,
    double * const p_mean_eps,
    double * const p_omega,
    double * const p_delta_psi
) {
    // TODO: explain
    const double mean_eps_arcsec = 84381.448 + jc * (-46.8150 + jc * (-0.00059 + jc * 0.001813));
    const double mean_eps = fmod(mean_eps_arcsec * ARCSEC_TO_RAD, 2 * PI);

    // Delaunay fundamental arguments in [deg]
    const double l_deg = 134.96298139 + jc * (1717915922.6330 + jc * (31.310 + jc * 0.064));
    const double l1_deg = 357.52772333 + jc * (129596581.2240 + jc * (-0.577 + jc * -0.012));
    const double f_deg = 93.27191028 + jc * (1739527263.1370 + jc * (-13.257 + jc * 0.011));
    const double d_deg = 297.85036306 + jc * (1602961601.3280 + jc * (-6.891 + jc * 0.019));
    const double omega_deg = 125.04452222 + jc * (-6962890.5390 + jc * (7.455 + jc * 0.008));

    // Convert from [deg] to [rad]
    const double l = deg_to_rad(fmod(l_deg, 360));
    const double l1 = deg_to_rad(fmod(l1_deg, 360));
    const double f = deg_to_rad(fmod(f_deg, 360));
    const double d = deg_to_rad(fmod(d_deg, 360));
    const double omega = deg_to_rad(fmod(omega_deg, 360));

    double delta_psi_p1mas = 0.0;
    double delta_eps_p1mas = 0.0;
    // Sum in reverse order to preserve floating point accuracy
    // TODO: see if this can be made size_t
    for (int i = LIB_GELOCUS_IAU80_NUTATION_TERMS - 1; i >= 0; --i)
    {
        const lib_gelocus_IAU80NutationCoeffSet * const coeff = &LIB_GELOCUS_IAU80_NUTATION_COEFFS[i];

        const double arg = (l * coeff->l) + (l1 * coeff->l1) + (f * coeff->f) + (d * coeff->d) + (omega * coeff->omg);

        delta_psi_p1mas += (coeff->sp + jc * coeff->spt) * sin(arg);
        delta_eps_p1mas += (coeff->ce + jc * coeff->cet) * cos(arg);
    }

    // Add in EOP corrections to GCRF
    delta_psi_p1mas += eop.dPsi;
    delta_eps_p1mas += eop.dEps;

    const double delta_psi = fmod(delta_psi_p1mas * P1MAS_TO_RAD, 2 * PI);
    const double delta_eps = fmod(delta_eps_p1mas * P1MAS_TO_RAD, 2 * PI);
    const double true_eps = mean_eps + delta_eps;

    const double sin_psi = sin(delta_psi);
    const double cos_psi = cos(delta_psi);
    const double sin_eps = sin(mean_eps);
    const double cos_eps = cos(mean_eps);
    const double sin_true_eps = sin(true_eps);
    const double cos_true_eps = cos(true_eps);

    if (N != NULL)
    {
        // Compute nutation matrix (TODO: explain which rotation)
        N->row1.x = cos_psi;
        N->row1.y = cos_true_eps * sin_psi;
        N->row1.z = sin_true_eps * sin_psi;
        N->row2.x = -cos_eps * sin_psi;
        N->row2.y = cos_true_eps * cos_eps * cos_psi + sin_true_eps * sin_eps;
        N->row2.z = sin_true_eps * cos_eps * cos_psi - cos_true_eps * sin_eps;
        N->row3.x = -sin_eps * sin_psi;
        N->row3.y = cos_true_eps * sin_eps * cos_psi - sin_true_eps * cos_eps;
        N->row3.z = sin_true_eps * sin_eps * cos_psi + cos_true_eps * cos_eps;
    }

    // Other optional outputs
    if (p_mean_eps != NULL)
    {
        *p_mean_eps = mean_eps;
    }
    if (p_omega != NULL)
    {
        *p_omega = omega;
    }
    if (p_delta_psi != NULL)
    {
        *p_delta_psi = delta_psi;
    }
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
