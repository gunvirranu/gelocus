<div class="title-block" style="text-align: center;" align="center">

# `🌐 lib_gelocus 🌐`

###### IAU-76 Earth-centered reference frame transformations in portable C99

[![Release][release_badge]][release]
[![License][license_badge]][license]
<!-- [![CI][ci_badge]][ci] -->

</div>

#### Note: WIP

`lib_gelocus` is a library for Earth-centered reference frames and transformations.
The base algorithms are written in strictly standards compliant C99 with no dependencies, with the aim of being simple, robust, portable, and usable for embedded or real-time uses.
Few assumptions are made and few restrictions are imposed, so it's relatively easy to adapt as needed for your use.

_State vectors_ represent time, position, and velocity in space, and are numerically realised relative to a reference frame.
A _reference frame_ is specified by an origin point (i.e. like the centre of the mass of the Earth) and an orientation relative to the universe that define its' coordinate axes.
There are some very commonly used frames used by the space community, such as those defined by the IAU-76/FK5 or IAU-2000 standards.
Given two reference frames, there are time-dependent linear transformations between them.
This library provides a simple interface for generating and applying these transformations.

The intention is that this can also easily serve as the underlying implementation for higher-level wrappers rather than us all redoing a bunch of work.
Towards this goal, a thin C++11 wrapper is provided to support a modern style of usage, with the coolest feature being templated compile-time tracking of reference frames (like SI units).
This means the type system itself will disallow accidentally mixing up state vectors in different frames, very cool.

## TODO

- CI pipeline
- Documentation:
  - installation
  - background
  - basic usage
  - examples
  - docstrings
- Benchmarking suite
- Impl new theories (IAU-2000 and beyond)

<!-- Links -->

<!-- Badges -->
[release]: https://github.com/gunvirranu/gelocus/releases "Release"
[release_badge]: https://img.shields.io/github/v/release/gunvirranu/gelocus?color=orange&logo=github&sort=semver "Release"
[license]: #license "License"
[license_badge]: https://img.shields.io/badge/license-MIT-blue.svg "License"
[ci]: https://github.com/gunvirranu/gelocus/actions/workflows/ci.yml "CI Workflow"
[ci_badge]: https://github.com/gunvirranu/gelocus/actions/workflows/ci.yml/badge.svg?branch=master "CI Workflow Badge"
