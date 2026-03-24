// VMD Plugin implemented as a Tcl dynamic library
// allows for loading of per-frame trajectory information from .toptraj files

#define PKG_NAME "toptraj"
#define VERSION "0.1"

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tcl8.6/tcl.h>
#include <tcl8.6/tclDecls.h>
#include <zconf.h>
#include <zlib.h>

/* Global state */
// sadly required, as the traces have to be reading information from somewhere

typedef struct {
    uint32_t frame_number;
    uint32_t n_atoms;
    // these are stored as tcl lists, in their string representation
    char *names;
    char *resnames;
    char *resids;
    char *types;
    char *charges;
    char *masses;
    char *bonds;
} TopTrajFrame;

typedef struct {
    // metadata of how it was made
    char *path;
    int molid;
    int remove_pbc_crossing;
    // data loaded into memory
    // number of frames loaded
    size_t n_frames;
    // current capacity of frames
    size_t cap_frames;
    TopTrajFrame **frames;
} TopTrajData;

static TopTrajData *load_toptraj(const char *path, int molid, int pbc) {
    /// Reads .toptraj file at path and loads it into a dynamically allocated
    /// object, which it returns.
    /// While doing so, prints a progress bar on STDOUT.
    /// Prints error message to STDERR and returns NULL if there is any error.
    TopTrajData *res = (TopTrajData *)calloc(sizeof(TopTrajData), 1);
    res->molid = molid;
    res->remove_pbc_crossing = pbc;
    // clone path, as we do not own it
    res->path = malloc(strlen(path) + 1);
    strcpy(res->path, path);
    // TODO
    return NULL;
}


static void free_frame(TopTrajFrame *frame) {
    free(frame->names);
    free(frame->resnames);
    free(frame->resids);
    free(frame->types);
    free(frame->charges);
    free(frame->masses);
    free(frame->bonds);
    free(frame);
}

static void free_toptraj(TopTrajData *data) {
    for (size_t i = 0; i < data->n_frames; i++) {
        free_frame(data->frames[i]);
    }
    if (data->frames != NULL) {
        free(data->frames);
    }
    if (data->path != NULL) {
        free(data->path);
    }
    free(data);
}

static char *on_frame_change(
    ClientData data, // TopTrajData *
    Tcl_Interp *interp,
    const char *name1, // "vmd_frame"
    const char *name2, // current frame as string
    int flags
) {
    // TODO
    return NULL;
}

static int load_toptraj_cmd(
    ClientData _data, // NULL
    Tcl_Interp *interp,
    int argc, char const *argv[]
) {
    if (argc < 2) {
        Tcl_SetResult(
            interp,
            "Usage: load_toptraj path_to_toptraj mol_id remove_pbc_crossing?",
            TCL_STATIC
        );
        return TCL_ERROR;
    }
    // which molecule ID is this trajectory for?
    int molid = 0;
    if (argc >= 3) {
        molid = atoi(argv[2]);
        if (strlen(argv[2]) == 0 || (molid == 0 && argv[2][0] != '0')) {
            Tcl_SetResult(
                interp,
                "molID should be a number",
                TCL_STATIC
            );
            return TCL_ERROR;
        }
    }
    // should we remove bonds that cross PBC?
    int pbc = 0;
    if (argc >= 4) {
        pbc = atoi(argv[3]);
    }
    TopTrajData *data = load_toptraj(argv[1], molid, pbc);
    Tcl_TraceVar(
        interp, "vmd_frame", TCL_TRACE_WRITES | TCL_GLOBAL_ONLY,
        on_frame_change, (ClientData)data
    );
    return TCL_OK;
}

static int unload_toptraj_cmd(
    ClientData _data, // NULL
    Tcl_Interp *interp,
    int argc, char const *argv[]
) {
    Tcl_SetResult(interp, "TODO: Unimplemented", TCL_STATIC);
    return TCL_ERROR;
}

int DLLEXPORT Toptraj_Init(Tcl_Interp *interp) {
    Tcl_PkgProvide(interp, PKG_NAME, VERSION);
    Tcl_CreateCommand(
        interp, "load_toptraj", load_toptraj_cmd, NULL, NULL
    );
    Tcl_CreateCommand(
        interp, "unload_toptraj", unload_toptraj_cmd, NULL, NULL
    );
    return TCL_OK;
}
