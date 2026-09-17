## VMD .toptraj plugin

This VMD plugin facilitates the visualization of dynamic bond networks, as well as it loads topology information
(atom names, types, charges, masses) into VMD in a dynamic manner. Dynamic in this context means that
different frames can have different values for it. This is achieved by setting a Tcl trace on the global
variable storing the VMD frame information to execute an efficient loader every time the frame is switched.

### Building

* Install dependencies (Ubuntu package names):
   - a C compiler and GNU Make (`build-essential`)
   - `pkg-config`
   - Tcl headers, same version as your VMD was compiled with (`tcl8.6-dev`)
   - 1g or ng-compat zlib headers (`zlib1g-dev`)
* Adjust Makefile if needed
* Run `make`
* Copy `toptraj.so` to a suitable location

### Usage

Load toptraj.so with the `load /path/to/toptraj.so` command in the Tcl console or script.
Loading it makes toptraj-specific Tcl commands available.

`load_toptraj path molID ?pbc?`

Params:

* path - path to .toptraj
* molID - molecule ID to use
* pbc - optional argument, if 1, it will ignore bonds that cross the PBC (only works for orthogonal PBC)
