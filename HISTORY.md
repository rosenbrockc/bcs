Revision History
======

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