# DPMJET-III

For current validation plots [click here](docs/figures/193/DPMJET-III-193-validation.pdf)

19.3.5
- rename any references to the code impy to chromo
- changed verbosity level to LPRi > 0 for logo
- changed some notation in Makefile

19.3.4:
- Different values for atomic mass of 12c depending on preprocessor flags
- Change to python interface of `pho_init`

19.3.3:

- Added preprcessor option for integration into CORSIKA 7
- Added preprcessor option for integration with impy that disables the internal random number generator
- Smaller fixes
- No change to physics

19.3.2:

- Build system adapted to different include schemes. Thx @luillo76 for PR
- New preprocessor flags: -DFLDOTINCL and -DFLINCINCL -PNUTINC
- No change to physics

19.3.1:

- Proprietary and nonfunctional FLUKA common blocks removed
- New method in `DT_GETPTN` to suppress short valence-(strange-sea) chains to suppress angular peaks in kaon distributions in low-energy AA interactions ([see this explanation](docs/validation.md))

19.2.1:

- No change to the model
- Added one example for photo hadronic interactions with PHOJET w/o using input cards

19.2.0:

- Synchronization of fragmentation parameters between this and the previous FLUKa version of DPMJET
- Minor but notable change to physics
- Explicit INT declaration in computed constant `NXZFBK` of (`FRBKCM`) for f2py

19.1.2:

- Minor changes of conditional compilations for FLUKA compatibility

19.1.1:

- Makefiles can build libraries with MinGW gfortran on Windows
- Executables don't work (yet) on Windows.
