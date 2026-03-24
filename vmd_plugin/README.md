VMD .toptraj plugin
===================

VMD plugin for loading .toptraj files.

Building
--------

1. Install dependencies dependencies (Ubuntu package names):
   - a C compiler and GNU Make (`build-essential`)
   - `pkg-config`
   - Tcl headers, same version as your VMD was compiled with (`tcl8.6-dev`)
    - Note: the Makefile is configured with system-wide installed Tcl 8.6 by default
   - "old" zlib headers (`zlib1g-dev`)
     - on some distros zlib-ng-compat could work, such as `zlib-ng-compat-devel` on Fedora
2. Run `make`

Usage
-----

1. Load toptraj.so with the `load /path/to/toptraj.so` command in the Tcl console or script.
2. Use the Tcl API to load the data (Data located in file out.toptraj, for molecule ID 0) `load_toptraj "out.toptraj" 0`.
