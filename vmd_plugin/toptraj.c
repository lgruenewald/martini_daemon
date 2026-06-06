// VMD Plugin implemented as a Tcl dynamic library
// allows for loading of per-frame trajectory information from .toptraj files

/*
TODO:
- selections other than all atoms
- deleting frames
- deleting traces
- everything here is prototype quality, fix that
    - use the Obj interface instead of string for tcl
    - non orthogonal pbc
*/

#define PKG_NAME "toptraj"
#define VERSION "0.1"

#include <assert.h>
#include <fcntl.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <tcl.h>
#include <tclDecls.h>
#include <zconf.h>
#include <zlib.h>

/* TYPES */

typedef struct bond {
    uint32_t bi;
    uint32_t bj;
} Bond;

typedef struct vec3 {
    float x;
    float y;
    float z;
} Vec3;

typedef struct chunk {
    size_t index;
    size_t len;
    char *content;
} Chunk;

typedef struct {
    bool is_err;
    char *error_msg;

    // metadata of how it was made
    char *path;
    int molid;
    int remove_pbc_crossing;

    // read from VMD trajectory
    size_t n_atoms;

    // data loaded into memory
    // mmap'd file
    int fd;
    size_t fsize;
    char *content;

    // per simulation data
    char *resnames;
    char *resids;

    // per frame data
    // number of frames loaded
    size_t n_frames;
    // current capacity of frames
    size_t cap_frames;
    // pointer to the start of frames on disk, into the mmap'd file
    char **frames;
} TopTrajData;

/// Decompress the buffer in, write decompressed to out.
static void decompress(
    TopTrajData *data,
    char *in, size_t in_len,
    char *out, size_t out_len,
    uint32_t crc32_expected
) {
    struct z_stream_s stream;
    stream.zalloc = Z_NULL;
    stream.zfree = Z_NULL;
    stream.opaque = Z_NULL;
    assert(in_len < INT32_MAX);
    assert(out_len < INT32_MAX);
    stream.avail_in = in_len;
    stream.avail_out = out_len;
    stream.next_in = (Bytef *)in;
    stream.next_out = (Bytef *)out;

    if (inflateInit(&stream) != Z_OK) {
        data->error_msg = "Couldn't initialize zlib stream.";
        data->is_err = true;
        return;
    }
    
    for (;;) {
        const int ret = inflate(&stream, Z_NO_FLUSH);
        switch (ret) {
            case Z_NEED_DICT:
            case Z_DATA_ERROR:
            case Z_STREAM_ERROR:
                data->error_msg = "Couldn't decompress. Incomplete or corrupt zlib compressed data.";
                data->is_err = true;
                goto finish;
            case Z_MEM_ERROR:
                data->error_msg = "Couldn't allocate memory.";
                data->is_err = true;
                goto finish;
            case Z_STREAM_END:
                goto finish;
            default:;
        }
    }

    finish:

    inflateEnd(&stream);

    if (data->is_err) {
        return;
    }

    uint32_t crc32_disk = crc32(0, (Bytef *)out, out_len);
    if (crc32_expected != crc32_disk) {
        data->error_msg = "CRC32 Mismatch. Toptraj file is corrupt.";
        data->is_err = true;
        return;
    }
}

static uint64_t parse_Q(Chunk *chunk) {
    uint64_t res;
    // assume little endian machine
    memcpy(&res, &chunk->content[chunk->index], 8);
    chunk->index += 8;
    return res;
}

static uint32_t parse_I(Chunk *chunk) {
    uint32_t res;
    // assume little endian machine
    memcpy(&res, &chunk->content[chunk->index], 4);
    chunk->index += 4;
    return res;
}

static float parse_f(Chunk *chunk) {
    float res;
    memcpy(&res, &chunk->content[chunk->index], 4);
    chunk->index += 4;
    return res;
}

static double parse_d(Chunk *chunk) {
    double res;
    memcpy(&res, &chunk->content[chunk->index], 8);
    chunk->index += 8;
    return res;
}

static char *read_s(Chunk *chunk) {
    const uint8_t size = chunk->content[chunk->index];
    const char *res = &chunk->content[chunk->index+1];
    chunk->index += 1 + size;
    char *s = malloc(size + 1);
    memcpy(s, res, size);
    s[size] = 0;
    return s;
}

#define READ_X(fname, f, type, fmt) \
static char *fname(Chunk *chunk, size_t n, long lim) { \
    type nums[n]; \
    size_t len = 0; \
    for (size_t i = 0; i < n; i++) { \
        nums[i] = f(chunk); \
        if (lim <= 0 || i < (size_t)lim) len += snprintf(NULL, 0, fmt, nums[i]); \
    } \
    char *res = malloc(len + 1); \
    size_t c = 0; \
    for (size_t i = 0; i < n; i++) { \
        if (lim <= 0 || i < (size_t)lim) c += sprintf(&res[c], fmt, nums[i]); \
    } \
    assert(c == len); \
    assert(res[len] == 0); \
    FREE \
    return res; \
}

#define FREE
READ_X(read_Is, parse_I, uint32_t, "%u ")
READ_X(read_fs, parse_f, float, "%.1f ")
#undef FREE
#define FREE \
for (size_t i = 0; i < n; i++) { \
    free(nums[i]); \
}
READ_X(read_ss, read_s, char *, "%s ")
#undef FREE

static bool sorted_list_insert(Bond *sorted_bonds, size_t *sorted_len, Bond b) {
    // find where to insert
    // highest index that is for sure valid
    size_t most = *sorted_len;
    // lowest index that is for sure valid
    size_t least = 0;

    while (most > least) {
        // rounding down intentional
        size_t guess_idx = (most + least) / 2;
        Bond guess = sorted_bonds[guess_idx];

        if (b.bi > guess.bi || (b.bi == guess.bi && b.bj > guess.bj)) {
            // bond is larger than guess
            // move least up
            least = guess_idx + 1;
        } else if (b.bi < guess.bi || (b.bi == guess.bi && b.bj < guess.bj)) {
            // bond is smaller than guess
            // move most down
            most = guess_idx;
            
        } else {
            // duplicate bond entry
            // forbidden by the file format, but we silently ignore it
            return false;
        }
        
    }

    // insert into the sorted list
    assert(most == least);
    if (least < *sorted_len) {
        memmove(&sorted_bonds[least+1], &sorted_bonds[least], (*sorted_len - least) * sizeof(Bond));
    }
    sorted_bonds[least] = b;
    (*sorted_len)++;
    return true;
    
}

static char *read_bonds(Chunk *chunk, size_t n_atoms, long lim, int pbc, Vec3 box, Vec3 *coords) {
    const size_t n_bonds = parse_Q(chunk);

    // sorted list of bonds
    size_t sorted_len = 0;
    Bond *sorted_bonds = calloc(2 * n_bonds, sizeof(struct bond));
    
    for (size_t i = 0; i < n_bonds; i++)
    {
        Bond b = { parse_I(chunk), parse_I(chunk) };
        if (b.bi == b.bj) {
            // self bonding
            // forbidden by the file format, but we silently ignore it
            continue;
        }
        if (b.bi >= n_atoms || b.bj >= n_atoms || (lim > 0 && (b.bi >= lim || b.bj >= lim))) {
            // one of the forming atoms is out of range
            // e.g. if solvent is removed from the trajectory
            continue;
        }
        if (pbc) {
            // skip if bond is longer than 1/2 along any pbc directions
            if (box.x / 2. < fabs(coords[b.bi].x - coords[b.bj].x)) {
                continue;
            }
            if (box.y / 2. < fabs(coords[b.bi].y - coords[b.bj].y)) {
                continue;
            }
            if (box.z / 2. < fabs(coords[b.bi].z - coords[b.bj].z)) {
                continue;
            }
        }
        // FIXME: good candidate for optimization, if it becomes a bottleneck
        sorted_list_insert(sorted_bonds, &sorted_len, b);
        Bond b2 = { b.bj, b.bi };
        sorted_list_insert(sorted_bonds, &sorted_len, b2);
    }

    size_t cap = 10000;
    size_t len = 0;
    char *res = malloc(cap);
    res[len] = '{';
    len++;

#define ACCOMODATE(x) \
    while (len+(x) > cap) { \
        cap *= 2; \
        res = realloc(res, cap); \
        assert(res != NULL); \
    }

    size_t bond_idx = 0;
    for (size_t i = 0; i < n_atoms; i++) {
        ACCOMODATE(1)
        res[len] = '{';
        len++;

        bool first_iter = true;
        while (bond_idx < sorted_len && sorted_bonds[bond_idx].bi == i) {
            if (!first_iter) {
                ACCOMODATE(1)
                res[len] = ' ';
                len++;
            }
            first_iter = false;

            Bond b = sorted_bonds[bond_idx];
            int extra = snprintf(NULL, 0, "%u", b.bj);
            ACCOMODATE(extra+1);
            snprintf(&res[len], extra+1, "%u", b.bj);
            len += extra;
            bond_idx++;
        }
        
        ACCOMODATE(2);
        res[len] = '}';
        len++;
        res[len] = ' ';
        len++;
    }

    ACCOMODATE(1);
    res[len] = '}';
    len++;

    res = realloc(res, len+1);
    assert(res != NULL);
    res[len] = 0;
    return res;
}

static Chunk *read_chunk(TopTrajData *data, Chunk *main) {
    /// Decompress a single chunk and return its contents.
    
    if (main->index >= main->len) {
        return NULL;
    }

    if (main->index + 24 >= main->len) {
        data->error_msg = "Corrupt .toptraj, leftover bytes found.";
        data->is_err = true;
        return NULL;
    }
    
    size_t len_comp = parse_Q(main);
    size_t len_decomp = parse_Q(main);
    char *content = malloc(len_decomp);
    char *comp = &main->content[main->index];
    main->index += len_comp;
    uint32_t crc32 = parse_I(main);
    
    decompress(
        data,
        comp, len_comp,
        content, len_decomp,
        crc32
    );

    Chunk *chunk = (Chunk *)calloc(1, sizeof(Chunk));
    chunk->content = content;
    chunk->len = len_decomp;
    chunk->index = 0;

    return chunk;
}

static void skip_chunk(TopTrajData *data, Chunk *main) {
    if (main->index >= main->len) {
        return;
    }

    if (main->index + 24 >= main->len) {
        data->error_msg = "Corrupt .toptraj, leftover bytes found.";
        data->is_err = true;
        return;
    }
    
    size_t len_comp = parse_Q(main);
    parse_Q(main); //len_decomp

    main->index += len_comp;
    parse_I(main); //crc32
}

static void free_chunk(Chunk *chunk) {
    free(chunk->content);
    free(chunk);
}

static void free_toptraj(TopTrajData *data) {
    if (data->frames != NULL) {
        free(data->frames);
    }
    if (data->path != NULL) {
        free(data->path);
    }
    if (data->content != NULL) {
        munmap(data->content, data->fsize);
    }
    if (data->resids != NULL) {
        free(data->resids);
    }
    if (data->resnames != NULL) {
        free(data->resnames);
    }
    free(data);
}

static TopTrajData *load_toptraj(Tcl_Interp *interp, const char *path, int molid, int pbc) {
    /// Reads .toptraj file at path and loads it into a dynamically allocated
    /// object, which it returns.
    /// While doing so, prints a progress bar on STDOUT.
    /// 
    /// Prints error message to STDERR and returns NULL if there is any error.
    TopTrajData *res = (TopTrajData *)calloc(1, sizeof(TopTrajData));
    res->is_err = false;

    res->molid = molid;
    res->remove_pbc_crossing = pbc;
    // clone path, as we do not own it
    res->path = strdup(path);

    res->fd = open(path, O_RDONLY);
    struct stat sb;
    fstat(res->fd, &sb);
    res->fsize = sb.st_size;
    res->content = mmap(
        NULL, res->fsize, PROT_READ, MAP_SHARED, res->fd, 0
    );
    if (res->content == MAP_FAILED) {
        res->is_err = true;
        res->error_msg = "Mapping the file into memory failed.";
        return res;
    }

    Chunk main;
    main.content = res->content;
    main.len = res->fsize;
    main.index = 0;

#define CHECK_ERR \
if (res->is_err) { \
    return res; \
}
    // 1.0 magic number
    if (strncmp(main.content, "\xc0TOPTR\x01\0", 8) == 0) {
        printf("TopTraj file version is 1.0\n");
        main.index += 8;

        // SIM NAME
        Chunk *header = read_chunk(res, &main);

        if (res->is_err) {
            return res;
        }

        // INITIAL MOLECULES
        char *name = read_s(header);
        printf("Simulation name: %s\n", name);
        const uint32_t n_initial_molecules = parse_I(header);
        printf("Number of initial molecule types: %u\n", n_initial_molecules);

        for (uint32_t i = 0; i < n_initial_molecules; i++)
        {
            char *molname = read_s(header);
            const uint32_t molcount = parse_I(header);
            const uint32_t atoms_per_mol = parse_I(header);
            printf("Mol name: %s\n", molname);
            printf("Mol count: %u\n", molcount);
            printf("Atoms per mol: %u\n", atoms_per_mol);
            free(molname);
        }

        int tcl_res;

        // RESNAME, RESID
        char *atomselect = NULL;
        char molid[64];
        snprintf(molid, 64, "%i", res->molid);
        tcl_res = Tcl_VarEval(interp, "atomselect ", molid, " all", NULL);
        if (tcl_res != TCL_OK) {
            printf("TOPTRAJ FATAL: failed at atomselect.\n");
        }
        atomselect = strdup(Tcl_GetStringResult(interp));
        tcl_res = Tcl_VarEval(interp, atomselect, " num", NULL);
        if (tcl_res != TCL_OK) {
            printf("TOPTRAJ FATAL: failed at num.\n");
        }
        long sel_atoms = atol(Tcl_GetStringResult(interp));

        uint32_t n_atoms = parse_I(header);
        res->resnames = read_ss(header, n_atoms, sel_atoms);
        res->resids = read_Is(header, n_atoms, sel_atoms);


        tcl_res = Tcl_VarEval(interp, atomselect, " set resname {", res->resnames, "}", NULL);
        if (tcl_res != TCL_OK) {
            printf("TOPTRAJ FATAL: failed at set resname.\n");
        }
        tcl_res = Tcl_VarEval(interp, atomselect, " set resid {", res->resids, "}", NULL);
        if (tcl_res != TCL_OK) {
            printf("TOPTRAJ FATAL: failed at set resid.\n");
        }
        tcl_res = Tcl_VarEval(interp, atomselect, " delete", NULL);
        if (tcl_res != TCL_OK) {
            printf("TOPTRAJ FATAL: failed at sel delete.\n");
        }

        free_chunk(header);
        free(atomselect);

        // SKIP OVER ALL FRAMES BUT SAVE OFFSETS
        res->cap_frames = 1000;
        res->frames = calloc(res->cap_frames, sizeof(char *));
        res->n_frames = 0;

        size_t frame_index = 0;
        printf("\n");

        while (main.index < main.len) {
            res->frames[frame_index] = &main.content[main.index];
            skip_chunk(res, &main);
            frame_index++;
            printf("Loaded frame %lu\r", frame_index);
            res->n_frames++;
            if (res->n_frames >= res->cap_frames)
            {
                res->cap_frames *= 2;
                res->frames = realloc(res->frames, sizeof(char *) * res->cap_frames);
                assert(res->frames != NULL);
            }
        }
        printf("\nRead %lu frames.\n", frame_index);

    } else {
        res->is_err = true;
        res->error_msg = "File is not a valid toptraj file, or the version is unknown.";
    }

    return res;
}



static char *on_frame_change(
    ClientData data, // TopTrajData *
    Tcl_Interp *interp,
    const char *name1, // "vmd_frame"
    const char *name2,
    int flags
) {
    // surpress unused warnings
    (void)flags;
    (void)name1;
    (void)name2;

    // all that needs freeing
    Chunk *chunk = NULL;
    char *names = NULL;
    char *types = NULL;
    char *charges = NULL;
    char *masses = NULL;
    char *bonds = NULL;
    char *sel = NULL;
    Vec3 *coords = NULL;
    
    TopTrajData *toptraj = (TopTrajData *)data;

    // GET CURRENT FRAME
    char molid[64];
    snprintf(molid, 64, "%i", toptraj->molid);

    // error handling
#define CHECK_RES(x) \
    if (res != TCL_OK) {\
        printf("TOPTRAJ FATAL: Trace failed at " x ": %s\n", Tcl_GetStringResult(interp)); \
        goto free; \
    }

#define CHECK_TR \
    if (toptraj->is_err) { \
        printf("%s\n", toptraj->error_msg); \
        toptraj->is_err = false; \
        toptraj->error_msg = NULL; \
        goto free; \
    }

    int res = Tcl_VarEval(interp, "molinfo ", molid, " get frame", NULL);
    CHECK_RES("molinfo get frame")
    int64_t c_frame = atol(Tcl_GetStringResult(interp));
    if (c_frame < 0 || (size_t)c_frame >= toptraj->n_frames) {
        return NULL;
    }

    char frame_str[64];
    snprintf(frame_str, 64, "%li", c_frame);

    // ATOMSELECT
    res = Tcl_VarEval(interp, "atomselect ", molid, " all frame ", frame_str, NULL);
    CHECK_RES("atomselect")

    sel = strdup(Tcl_GetStringResult(interp));

    // NUMBER OF ATOMS IN SEL
    res = Tcl_VarEval(interp, sel, " num", NULL);
    CHECK_RES("failed at num")
    long sel_atoms = atol(Tcl_GetStringResult(interp));

    // DECOMPRESS FRAME DATA
    Chunk main;
    main.content = toptraj->frames[c_frame];
    main.index = 0;
    main.len = toptraj->fsize - (main.content - toptraj->content);

    chunk = read_chunk(toptraj, &main);
    CHECK_TR

    uint32_t frame_number = parse_I(chunk);
    (void)frame_number;
    CHECK_TR
    uint32_t n_atoms = parse_I(chunk);
    CHECK_TR
    size_t sim_step = parse_Q(chunk);
    (void)sim_step;
    CHECK_TR
    double sim_time_ns = parse_d(chunk);
    (void)sim_time_ns;
    CHECK_TR

    names = read_ss(chunk, n_atoms, sel_atoms);
    CHECK_TR
    types = read_ss(chunk, n_atoms, sel_atoms);
    CHECK_TR
    charges = read_fs(chunk, n_atoms, sel_atoms);
    CHECK_TR
    masses = read_fs(chunk, n_atoms, sel_atoms);
    CHECK_TR

    // REMOVE PBC BONDS IF NEEDED
    Vec3 box = {0., 0., 0.};
    // if ignoring bonds that cross pbc, read out the coords
    if (toptraj->remove_pbc_crossing) {
        res = Tcl_VarEval(interp, sel, " get x", NULL);
        if (res != TCL_OK) {
            printf("TOPTRAJ FATAL: failed at get x: %s\n", Tcl_GetStringResult(interp));
            goto free;
        }
        char *xs = strdup(Tcl_GetStringResult(interp));
        res = Tcl_VarEval(interp, sel, " get y", NULL);
        if (res != TCL_OK) {
            printf("TOPTRAJ FATAL: failed at get y: %s\n", Tcl_GetStringResult(interp));
            free(xs);
            goto free;
        }
        char *ys = strdup(Tcl_GetStringResult(interp));
        res = Tcl_VarEval(interp, sel, " get z", NULL);
        if (res != TCL_OK) {
            printf("TOPTRAJ FATAL: failed at get z: %s\n", Tcl_GetStringResult(interp));
            free(xs);
            free(ys);
            goto free;
        }
        char *zs = strdup(Tcl_GetStringResult(interp));

        coords = calloc(n_atoms, sizeof(Vec3));
        char *xp = xs;
        char *yp = ys;
        char *zp = zs;
        for (size_t j = 0; j < n_atoms; j++) {
            // strtof autoskips whitespace
            Vec3 v = {
                strtof(xp, &xp),
                strtof(yp, &yp),
                strtof(zp, &zp)
            };
            coords[j] = v;
        }

        res = Tcl_VarEval(interp, "molinfo ", molid, " get {a b c}", NULL);
        if (res != TCL_OK) {
            printf("TOPTRAJ FATAL: failed at get pbc: %s\n", Tcl_GetStringResult(interp));
            free(xs);
            free(ys);
            free(zs);
            goto free;
        }
        char *bs = strdup(Tcl_GetStringResult(interp));
        char *bp = bs;
        box.x = strtof(bp, &bp);
        box.y = strtof(bp, &bp);
        box.z = strtof(bp, &bp);
        free(bs);
        free(xs);
        free(ys);
        free(zs);
    }
    // READ BONDS
            
    bonds = read_bonds(
        chunk, n_atoms, sel_atoms,
        toptraj->remove_pbc_crossing, box, coords
    );

    // SET ATOM PROPERTIES

    res = Tcl_VarEval(interp, sel, " set name {", names, "}", NULL);
    CHECK_RES("set name")
    res = Tcl_VarEval(interp, sel, " set type {", types, "}", NULL);
    CHECK_RES("set type")
    res = Tcl_VarEval(interp, sel, " set charge {", charges, "}", NULL);
    CHECK_RES("set charge")
    res = Tcl_VarEval(interp, sel, " set mass {", masses, "}", NULL);
    CHECK_RES("set mass")

    // SETBONDS
    res = Tcl_VarEval(interp, sel, " setbonds ", bonds, NULL);
    CHECK_RES("setbonds")

    res = Tcl_VarEval(interp, sel, " delete", NULL);
    CHECK_RES("delete")

free:
    if (chunk != NULL) free_chunk(chunk);
    if (names != NULL) free(names);
    if (types != NULL) free(types);
    if (charges != NULL) free(charges);
    if (masses != NULL) free(masses);
    if (bonds != NULL) free(bonds);
    if (sel != NULL) free(sel);
    if (coords != NULL) free(coords);

    return NULL;
}

static int load_toptraj_cmd(
    ClientData _data, // NULL
    Tcl_Interp *interp,
    int argc, char const *argv[]
) {
    (void)_data; // supress unused warning
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
    TopTrajData *data = load_toptraj(interp, argv[1], molid, pbc);
    if (data->is_err) {
        Tcl_SetResult(
            interp,
            data->error_msg,
            TCL_STATIC
        );
        free_toptraj(data);
        return TCL_ERROR;
    }
    Tcl_TraceVar(
        interp, "vmd_frame", TCL_TRACE_WRITES | TCL_GLOBAL_ONLY,
        on_frame_change, (ClientData)data
    );
    return TCL_OK;
}

int DLLEXPORT Toptraj_Init(Tcl_Interp *interp) {
    Tcl_PkgProvide(interp, PKG_NAME, VERSION);
    Tcl_CreateCommand(
        interp, "load_toptraj", load_toptraj_cmd, NULL, NULL
    );
    return TCL_OK;
}
