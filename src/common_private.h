/*
 * gelocus -- Earth-centered reference frame transformations library
 *
 * "Scientific code that shouldn't feel like it was written by scientists!"
 *
 * Licensed under the MIT License <http://opensource.org/licenses/MIT>.
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2024 Gunvir Singh Ranu
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef GELOCUS_COMMON_PRIVATE_H
#define GELOCUS_COMMON_PRIVATE_H

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
