## DPMJET-III

19.1.6:
- Removed old python lib build from Makefile
- Added new `ISWmdl(6) = 4` flag forcing Pythia to accept decay settings defined in mdcy 

19.1.5:
- Renamed impy references to chromo due to name change of the package 

19.1.4:

- No change to physics or code. Changed python interface definition of `pho_init`

19.1.3:

- No change to physics or code. Only conditional compilation added for `impy`

19.1.2:

- Minor changes of conditional compilations for FLUKA compatibility

19.1.1:

- Makefiles can build libraries with MinGW gfortran on Windows
- Executables don't work (yet) on Windows.
