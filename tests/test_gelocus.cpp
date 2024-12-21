#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

extern "C" {
#  include <gelocus/gelocus.h>
};

// Save some typing
#define jd_to_jc                lib_gelocus_jd_to_jc
typedef lib_gelocus_Vec3        Vec3;
typedef lib_gelocus_Matrix3     Matrix3x3;

TEST_CASE("test_constants")
{
    // Check J2000 epochs in julian date and julian century are the same
    CHECK(jd_to_jc(LIB_GELOCUS_EPOCH_J2000_JD) == LIB_GELOCUS_EPOCH_J2000_JC);

    // Check 1 day delta in julian date is same as in julian century
    const double jc_next_day = lib_gelocus_jd_to_jc(LIB_GELOCUS_EPOCH_J2000_JD + LIB_GELOCUS_DELTA_JD_PER_DAY);
    CHECK(jc_next_day == (LIB_GELOCUS_EPOCH_J2000_JC + LIB_GELOCUS_DELTA_JC_PER_DAY));
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

}

TEST_CASE("test_jd_to_jc")
{

}

TEST_CASE("test_apply_transformation")
{

}

TEST_CASE("test_transform")
{

}

TEST_CASE("test_compute_transformation_matrix")
{

}
