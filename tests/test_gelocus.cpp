/**
 * SPDX-License-Identifier: MIT
 * Copyright (c) 2024 Gunvir Singh Ranu
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

extern "C" {
#  include <lib_gelocus/gelocus.h>
#  include "../src/common_private.h"
};

using doctest::Approx;

// Save some typing
typedef lib_gelocus_Vec3        Vec3;
typedef lib_gelocus_Matrix3     Matrix3;

/// Useful macro to check `Vec3` equality
#define CHECK_VEC(lhs, rhs, eps) \
    do { \
        CHECK(lhs.x == Approx(rhs.x).epsilon(eps)); \
        CHECK(lhs.y == Approx(rhs.y).epsilon(eps)); \
        CHECK(lhs.z == Approx(rhs.z).epsilon(eps)); \
    } while (false)

TEST_CASE("test_constants")
{
    // Check pi is behaving like pi
    CHECK(PI == Approx(3).epsilon(0.05));
    CHECK(RAD_TO_DEG * DEG_TO_RAD == 1);

    // Check arcmin and arcsec conversions make sense
    CHECK(DEG_TO_ARCMIN * ARCMIN_TO_ARCSEC == DEG_TO_ARCSEC);
    CHECK(DEG_TO_RAD / DEG_TO_ARCMIN == Approx(ARCMIN_TO_RAD).epsilon(1e-16));
    CHECK(DEG_TO_RAD / DEG_TO_ARCSEC == Approx(ARCSEC_TO_RAD).epsilon(1e-16));
    CHECK(1e-3 * ARCSEC_TO_RAD == Approx(MAS_TO_RAD).epsilon(1e-16));
    CHECK(1e-4 * ARCSEC_TO_RAD == Approx(P1MAS_TO_RAD).epsilon(1e-16));

    // Check J2000 epochs in julian date and julian century are the same
    CHECK(lib_gelocus_jd_to_jc(LIB_GELOCUS_EPOCH_J2000_JD) == LIB_GELOCUS_EPOCH_J2000_JC);

    // Check 1 day delta in julian date is same as in julian century
    const double jc_next_day = lib_gelocus_jd_to_jc(LIB_GELOCUS_EPOCH_J2000_JD + LIB_GELOCUS_DELTA_JD_PER_DAY);
    CHECK(jc_next_day == (LIB_GELOCUS_EPOCH_J2000_JC + LIB_GELOCUS_DELTA_JC_PER_DAY));

    // Check J1900 epoch is roughly -1 centuries back from J2000
    CHECK(lib_gelocus_jd_to_jc(LIB_GELOCUS_EPOCH_J1900_JD) == Approx(-1).epsilon(1e-4));

    // Check GPS epoch (~1980) roughly matches fractional centuries
    CHECK(lib_gelocus_jd_to_jc(LIB_GELOCUS_EPOCH_GPS_JD) == Approx(-0.2).epsilon(1e-3));

    // Check day and JC per day are inverses
    CHECK((1 / LIB_GELOCUS_DELTA_JC_PER_DAY) == Approx(LIB_GELOCUS_DELTA_DAY_PER_JC).epsilon(1e-20));

    // Check identity matrix elements
    CHECK(LIB_GELOCUS_MATRIX_IDENTITY.row1.x == 1);
    CHECK(LIB_GELOCUS_MATRIX_IDENTITY.row1.y == 0);
    CHECK(LIB_GELOCUS_MATRIX_IDENTITY.row1.z == 0);
    CHECK(LIB_GELOCUS_MATRIX_IDENTITY.row2.x == 0);
    CHECK(LIB_GELOCUS_MATRIX_IDENTITY.row2.y == 1);
    CHECK(LIB_GELOCUS_MATRIX_IDENTITY.row2.z == 0);
    CHECK(LIB_GELOCUS_MATRIX_IDENTITY.row3.x == 0);
    CHECK(LIB_GELOCUS_MATRIX_IDENTITY.row3.y == 0);
    CHECK(LIB_GELOCUS_MATRIX_IDENTITY.row3.z == 1);
}

TEST_CASE("test_vec_norm")
{
    // Check default is 0
    Vec3 x = { 0 };
    CHECK(lib_gelocus_vec_norm(x) == 0);

    // Check basic pythagorean quadruple
    x = { 4, 13, 16 };
    CHECK(lib_gelocus_vec_norm(x) == 21);

    // Check norm is unchanged by negatives
    x = { -4, 13, -16 };
    CHECK(lib_gelocus_vec_norm(x) == 21);

    // Check all negatives and some floating point
    x = { -2, -3, -4 };
    CHECK(lib_gelocus_vec_norm(x) == Approx(5.3851648071).epsilon(1e-10));
}

TEST_CASE("test_dot_product")
{
    Vec3 a = { 0 };
    Vec3 b = { 0 };

    // Check default is 0
    CHECK(lib_gelocus_dot_product(a, b) == 0);

    // Check nominal dot product
    a = { 1, 2, 3 };
    b = { 4, 5, 6 };
    CHECK(lib_gelocus_dot_product(a, b) == 32);

    // Check commutative
    CHECK(lib_gelocus_dot_product(b, a) == lib_gelocus_dot_product(a, b));

    // Check negatives apply
    a = { 1, -2, -4};
    b = { 4, 5, -5 };
    CHECK(lib_gelocus_dot_product(a, b) == 14);
}

TEST_CASE("test_multiply_matrix")
{
    Matrix3 A = { 0 };
       Vec3 x = { 0 };
       Vec3 b = { 0 };
       Vec3 truth = { 0 };

    // Check default values are all 0
    b = lib_gelocus_multiply_matrix(A, x);
    CHECK_VEC(b, truth, 1e-30);

    // Check identity matrix does nothing
    A = LIB_GELOCUS_MATRIX_IDENTITY;
    x = { 7, -3, 12 };
    b = lib_gelocus_multiply_matrix(A, x);
    truth = x;
    CHECK_VEC(b, truth, 1e-30);

    // Check axes shuffling
    A = { .row1 = { 0,  0, -1 },
          .row2 = { 1,  0,  0 },
          .row3 = { 0, -1,  0 } };
    x = { -5, 9, 4 };
    b = lib_gelocus_multiply_matrix(A, x);
    truth = { -4, -5, -9 };
    CHECK_VEC(b, truth, 1e-30);

    // Test some random rotation matrix
    A = { .row1 = {  0.87559502, -0.2385524 ,  0.42003109 },
          .row2 = {  0.29597008,  0.95215193, -0.07621294 },
          .row3 = { -0.38175263,  0.19104831,  0.90430386 } };
    x = { 197.09117549, 51.11557365, 230.36351425 };
    b = lib_gelocus_multiply_matrix(A, x);
    truth = { 257.13814674, 89.44620389, 142.84408326 };  // Via faith in NumPy
    CHECK_VEC(b, truth, 1e-7);
}

TEST_CASE("test_jd_to_jc")
{
    // Vallado Example 3-5
    CHECK(lib_gelocus_jd_to_jc(2448855.009722) == Approx(-0.073647919).epsilon(1e-8));

    // Vallado Example 3-7
    CHECK(lib_gelocus_jd_to_jc(2453140.19727065) == Approx(0.043674121031).epsilon(1e-12));

    // Vallado Example 3-7
    CHECK(lib_gelocus_jd_to_jc(2453140.196522415) == Approx(0.043671100545).epsilon(1e-4));
}

TEST_CASE("test_jd_frac_to_jc")
{
    // Vallado Example 3-5
    CHECK(lib_gelocus_jd_frac_to_jc(2448855, 0.009722) == Approx(-0.073647919).epsilon(1e-9));

    // Vallado Example 3-7
    CHECK(lib_gelocus_jd_frac_to_jc(2453140, 0.19727065) == Approx(0.043674121031).epsilon(1e-12));

    // Vallado Example 3-7
    CHECK(lib_gelocus_jd_frac_to_jc(2453140, 0.196522415) == Approx(0.043671100545).epsilon(1e-5));
}

TEST_CASE("test_apply_transformation")
{
    lib_gelocus_Transformation trans = {
        .jc = LIB_GELOCUS_EPOCH_J2000_JC,
        .from = LIB_GELOCUS_FRAME_J2000,
        .to   = LIB_GELOCUS_FRAME_ECEF,
        .eop = {
            .dPsi = 0.01,
        },
        .matrix = LIB_GELOCUS_MATRIX_IDENTITY,
    };
    lib_gelocus_StateVector in = {
        .jc = LIB_GELOCUS_EPOCH_J2000_JC,
        .pos = { 7, 8, 9 },
        .vel = { 0, 0, 0 },
        .eop = { 0 },
        .frame = LIB_GELOCUS_FRAME_J2000,
    };
    lib_gelocus_StateVector out = { 0 };
    Vec3 vec = { 0 };

    // Verify no work if NULL pointer
    CHECK_FALSE(lib_gelocus_apply_transformation(trans, in, NULL));

    // Verify no work if from frame doesn't match input vector
    in.frame = LIB_GELOCUS_FRAME_ECEF;
    CHECK_FALSE(lib_gelocus_apply_transformation(trans, in, &out));
    // Verify `out` EOP data is untouched
    CHECK(out.eop.dPsi == 0);
    // Verify `out` vector is untouched
    CHECK_VEC(out.pos, vec, 1e-20);
    // Revert
    in.frame = LIB_GELOCUS_FRAME_J2000;

    // Verify no work if timestamp mismatch
    in.jc += (100 * LIB_GELOCUS_DELTA_JC_PER_DAY);
    CHECK_FALSE(lib_gelocus_apply_transformation(trans, in, &out));
    // Verify `out` EOP data is untouched
    CHECK(out.eop.dPsi == 0);
    // Verify `out` vector is untouched
    CHECK_VEC(out.pos, vec, 1e-20);
    // Revert
    in.jc = LIB_GELOCUS_EPOCH_J2000_JC;

    // Verify transformation is applied if conditions correct
    CHECK(lib_gelocus_apply_transformation(trans, in, &out));
    // Verfy output state vector timestamp is copied
    CHECK(out.jc == LIB_GELOCUS_EPOCH_J2000_JC);
    // Verify output frame makes sense
    CHECK(out.frame == trans.to);
    // Verify EOP data from transformation is copied
    CHECK(out.eop.dPsi == trans.eop.dPsi);
    // Verify output vec is computed (matrix is identity)
    vec = in.pos;
    CHECK_VEC(out.pos, vec, 1e-20);
}

TEST_CASE("test_transform")
{
    lib_gelocus_StateVector in = {
        .jc = LIB_GELOCUS_EPOCH_J2000_JC,
        .pos = { 300, 200, 400 },
        .vel = { 0, 0, 0 },
        .eop = { 0 },
        .frame = LIB_GELOCUS_FRAME_J2000,
    };

    // Verify output frame and timestamp make sense
    const lib_gelocus_StateVector out = lib_gelocus_transform(in, LIB_GELOCUS_FRAME_ECEF);
    CHECK(out.frame == LIB_GELOCUS_FRAME_ECEF);
    CHECK(out.jc == in.jc);

    // TODO: test actual transformation values
}

TEST_CASE("test_compute_transformation_matrix" * doctest::skip())
{

}
