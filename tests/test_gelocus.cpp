#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

extern "C" {
#  include <gelocus/gelocus.h>
};

TEST_CASE("test_constants")
{
    CHECK(lib_gelocus_jd_to_jc(LIB_GELOCUS_EPOCH_J2000_JD) == LIB_GELOCUS_EPOCH_J2000_JC);

    const double jc_next_day = lib_gelocus_jd_to_jc(LIB_GELOCUS_EPOCH_J2000_JD + LIB_GELOCUS_DELTA_JD_PER_DAY);
    CHECK(jc_next_day == (LIB_GELOCUS_EPOCH_J2000_JC + LIB_GELOCUS_DELTA_JC_PER_DAY));
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
    CHECK(3 == (4 + 1));
}
