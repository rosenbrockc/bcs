Revision History
======

Revision 1.0.6
------

- Added python wrapping support for _some_ of the functions. Those which don't have any `allocatable, intent(out)` parameters can be compiled into a python module my calling `make F90=gfortran` in the new `./wrap/` directory. It generates a shared library `_bcs.so` and a python module `bcs.py`. In order to run correctly, it needs [f90wrap](https://github.com/jameskermode/f90wrap) installed.
- Changed the function `file_exists` to use `inquire` instead of attempting an open since that creates the missing file on some compilers.
- Added a primitive control for the oscillation problem where the solution jumps between two `l0` norms.

Revision 1.0.5
------

- Worked around a bug in `gfortran` where assignment of a variable's value inside of an associate caused spurious errors with subtraction not being performed correctly. We have yet post the bug to `gfortran`'s bug list.

Revision 1.0.4
------

- Fixed a segfault with allocating optional `seed` parameters for subroutines that use random numbers.
- Fixed a segfault with calculating RMS and absolute error for the holdout set.

Revision 1.0.3
------

- Fixed the unit tests for the `do_bcs()`.
- Added the write-up explaining the mathematics, though it still needs a table updated.

Revision 1.0.2
------

- Updated the README file with citation request and information for contributing.

Revision 1.0.1
------

- Added the data files for all the unit tests.

Revision 1.0
------

- Finished documenting and coding all the modules.
- Added unit tests for the entire `bcs.f90` module. There appears to be a bug with the output from `do_bcs()` where only the first vector in each fit is being written to file. The numbers seem reasonable otherwise. I still need to confirm the output of that routine with Mathematica.