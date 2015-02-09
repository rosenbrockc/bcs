Revision History
======

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