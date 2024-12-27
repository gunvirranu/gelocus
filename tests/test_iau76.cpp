#include <doctest/doctest.h>

extern "C" {
#  include <gelocus/gelocus.h>
#  include <gelocus/iau76.h>
};

TEST_CASE("test_gmst")
{
    CHECK(lib_gelocus_gmst(LIB_GELOCUS_EPOCH_J2000_JC) > 0);
}

TEST_CASE("test_iau76_precession" * doctest::skip())
{

}
