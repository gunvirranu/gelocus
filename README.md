# gelocus

## Note: WIP

`gelocus` (from *ge locus*) is a C and C++ library for Earth-centered reference frame transformations.
_State vectors_ represent time, position, and velocity in space, and are numerically realised relative to a reference frame.
A _reference frame_ is specified by an origin point (i.e. like the centre of the mass of the Earth) and an orientation relative to the universe that define its' coordinate axes.
Given two reference frames, there are time-dependent linear transformations between them.
This library provides a simple interface for generating and applying these transformations.

The `gelocus` base lib is written in strictly standards compliant C99 with no dependencies, with the aim of being simple, robust, portable, and usable in real-time embedded contexts.
A thin C++11 wrapper is provided to support a modern style of usage, with the coolest feature being templated comile time tracking of reference frames (like SI units).
This means the type system itself will disallow accidentally mixing up state vectors in different frames, very cool.

## TODO:

- Code comments
- Proper documentation, installation, usage, examples, general details
- Tests
- Benchmarks
- Look into implementing newer theories (IAU-2000 and beyond)

<!-- Links -->
