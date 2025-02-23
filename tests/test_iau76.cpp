/**
 * SPDX-License-Identifier: MIT
 * Copyright (c) 2024 Gunvir Singh Ranu
 */

#include <doctest/doctest.h>

extern "C" {
#  include <lib_gelocus/gelocus.h>
#  include <lib_gelocus/iau76.h>
};

TEST_CASE("test_gmst")
{
    CHECK(lib_gelocus_gmst(LIB_GELOCUS_EPOCH_J2000_JC) > 0);
}

TEST_CASE("test_iau76_precession" * doctest::skip())
{

}
