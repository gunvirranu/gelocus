#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

extern "C" {
#  include <gelocus/gelocus.h>
};

TEST_CASE("test_foo")
{
    CHECK(1 == 1);
    CHECK(3 == 4);
}