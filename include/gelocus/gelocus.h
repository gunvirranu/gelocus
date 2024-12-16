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

#ifndef GELOCUS_GELOCUS_H
#define GELOCUS_GELOCUS_H

#include <stdbool.h>

// TODO: explain J2000 epoch meaning
const double LIB_GELOCUS_EPOCH_J2000_JD = 2451545.0;      ///< [julian date]    Time stamp of J2000 epoch
const double LIB_GELOCUS_EPOCH_J2000_JC = 0.0;            ///< [julian century] Time stamp of J2000 epoch
const double LIB_GELOCUS_DELTA_JD_PER_DAY = 1.0;          ///< [julian date]    Delta per 24 hour day
const double LIB_GELOCUS_DELTA_JC_PER_DAY = 1.0 / 36525;  ///< [julian century] Delta per 24 hour day

/// A 3-dimensional vector
typedef struct {
    double x;  ///< Index [0]
    double y;  ///< Index [1]
    double z;  ///< Index [2]
} lib_gelocus_Vec3;

/// A 3×3 matrix for linear transformations
///
/// If your larger project doesn't already vector/matrix types, feel free to just use
/// [lib_gelocus_StateVector] as your native type. However, if the larger project you're
/// integrating into uses some other vector/matrix types, I can't think of a great way of
/// supporting using those here without a bunch of hacks. My first suggestion is that you
/// just shuttle these types back and forth as needed at the interface layer, wherever
/// you call/wrap `gelocus` functions. If you truly need this library to use your native
/// types (e.g. you *really* need this to be vectorized), I am afraid you will need to
/// dig into this source and modify as needed to use your types. It should be relatively
/// straightforward. Let me know if you think of a cool way to support both.
typedef struct {
    lib_gelocus_Vec3 col1;  ///< Column 1 of 3×3 matrix, indices [i, 0]
    lib_gelocus_Vec3 col2;  ///< Column 2 of 3×3 matrix, indices [i, 1]
    lib_gelocus_Vec3 col3;  ///< Column 3 of 3×3 matrix, indices [i, 2]
} lib_gelocus_Matrix3;

/// Reference Frames
///
/// A reference frame defines coordinates on the affine space that is actual space.
/// It comprises an origin point and a 3D orientation defining the axes.
/// The origin for all frames considered in this library, which happens to be
/// prefixed with *ge*, is the centre of mass (CoM) of the Earth. This mean the
/// zero vector [0, 0, 0] in all frames is Earth's CoM. The x, y, and z axes are
/// the ones that change across the below frames.
///
/// The primary effort of this library is defining the numerical transformations
/// between the following frames, allowing simple conversion and computation.
/// These essentially boil down to just linear transformations, but the tricky
/// part is relations between depend on time. For example, the relation between
/// the ECI J2000 frame (which is meant to be inertial and non-moving relative
/// to the external universe) and the ECEF frame (which is fixed to Earth's surface)
/// changes as the Earth rotates.
typedef enum {
    LIB_GELOCUS_FRAME_J200,     ///< ECU: Earth-Centered Inertial fixed to the J2000 epoch
    LIB_GELOCUS_FRAME_MOD,      ///< MOD: Mean of Date (precession applied)
    LIB_GELOCUS_FRAME_TOD,      ///< TOD: True of Data (nutation applied)
    LIB_GELOCUS_FRAME_TEME,     ///< TEME: True Equator Mean Equinox, used only for SGP4
    LIB_GELOCUS_FRAME_PEF,      ///< PEF: Pseudo-Earth fixed (TODO: explain)
    LIB_GELOCUS_FRAME_ECEF      ///< ECEF: Earth-Centered Earth-Fixed
} lib_gelocus_Frame;

/// Earth Orientation Parameters
///
/// These are measured statistical corrections to the transformation models.
/// They are published by external sources and improve accuracy of conversions.
/// You may leave them as default (all 0) if you don't wish to provide these updates.
typedef struct {
    double xp;      ///< [rad]
    double yp;      ///< [rad]
    double dPsi;
    double dEps;
    double dX;
    double dY;      ///< Not used for IAU-76/FK5 transformations
} lib_gelocus_EOPData;

/// Space State Vector
///
/// This represents a dynamical/kinematic state of an object in space relative
/// to a reference frame. It contains a time, position, velocity, any extra EOP
/// parameters, and of course, the reference frame itself.
typedef struct {
    double jc;                  ///< [julian century] Time stamp
    lib_gelocus_Vec3    pos;    ///< [km] Position vector
    lib_gelocus_Vec3    vel;    ///< [km/s] Velocity vector
    lib_gelocus_Frame   frame;  ///< Reference frame which coordinates are relative to
    lib_gelocus_EOPData eop;    ///< EOP correction parameters at this time (if any)
} lib_gelocus_StateVector;

/// Numerical transformation between two reference frames
///
/// Valid per specific choice of origin and destination frame, specific time-stamp, and
/// and any given EOP parameters. If you're just looking to convert state vectors, prefer
/// to use [lib_gelocus_transform]. Getting an explicit transformation is mainly useful
/// if you want to save it for re-use rather than re-computing many times. In this case,
/// you may specify the time `jc`, the `from` and `to` frames, provide any `eop` data,
/// and then use [lib_gelocus_compute_transformation_matrix] to compute the `matrix`.
/// It can then be later applied to a state vector with [lib_gelocus_apply_transformation].
typedef struct {
    double jc;
    lib_gelocus_Frame   from;
    lib_gelocus_Frame   to;
    lib_gelocus_EOPData eop;
    lib_gelocus_Matrix3 matrix;
} lib_gelocus_Transformation;

/// Multiply a 3×3 matrix by a 3×1 vector, returning a 3×1 vector
lib_gelocus_Vec3 lib_gelocus_multiply_matrix(lib_gelocus_Matrix3 mat, lib_gelocus_Vec3 vec);

/// Convert a [julian date] to a [julian century].
///
/// TODO: consider sanity check on conversion?
double lib_gelocus_jd_to_jc(double jd);

/// Apply a transformation (between two frames) to a state vector
///
/// Returns `true` if the transformation is valid and achieved.
/// Returns `false` if invalid (i.e. origin frame of given [lib_gelocus_Transformation]
/// does not match frame of [lib_gelocus_StateVector]). This can happen if you mixed up
/// frames between generating the transformation and applying to a state vector.
bool lib_gelocus_apply_transformation(
    lib_gelocus_Transformation trans,
    lib_gelocus_StateVector sv,
    lib_gelocus_StateVector * out
);

/// Transform a state vector in one frame to another
lib_gelocus_StateVector lib_gelocus_transform(lib_gelocus_StateVector sv, lib_gelocus_Frame to);

/// Generate the 3 x 3 transformation between two frames. It also relies on the time
/// and the current EOP corrections. The 9-length matrix array contains the row-major
/// indices of the transformation. You shouldn't need to use this method unless you
/// want to peek at the internal matrix yourself; idk why, your call.
///
/// Assumes *all* fields (`jc`, `from`, `to`, `eop`) of [lib_gelocus_Transformation] are
/// populated. This only computes/fills out the `matrix` field, previous values unused.
void lib_gelocus_compute_transformation_matrix(lib_gelocus_Transformation * trans);

#endif  // GELOCUS_GELOCUS_H
