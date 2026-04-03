// VMD Plugin implemented as a Tcl dynamic library
// allows for loading of per-frame trajectory information from .toptraj files

/*
TODO:
- selections other than all atoms
- deleting frames
*/

#define PKG_NAME "toptraj"
#define VERSION "0.1"

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <tclDecls.h>
#include <zconf.h>
#include <zlib.h>

/* Streaming compressed file reader */
typedef struct
{
    // error handling
    bool is_err;
    char *error_msg; // always static string

    FILE *file;

    char *header;
    struct z_stream_s stream;

    size_t chunk_size;
    Bytef *chunk_in;
    Bytef *chunk_out;

    // within chunk_out, where are we at
    size_t read_out;

    bool stream_over;

    // current value of crc
    uLong crc;
} CompressedReader;


/// reads as much of the file into chunk as possible into the chunk
/// if stream.avail_in > 0, copies the leftovers first
/// sets stream.next_in and stream.avail_in
static void read_file_into_chunk(CompressedReader *reader) {
    if (reader->stream.avail_in > 0) {
        memmove(
            reader->chunk_in,
            reader->stream.next_in,
            reader->stream.avail_in
        );
    }
    const size_t read_from_file = reader->chunk_size - reader->stream.avail_in;
    const size_t read = fread(&reader->chunk_in[reader->stream.avail_in], 1, read_from_file, reader->file);
    if (read == 0) {
        reader->stream_over = true;
    }
    reader->stream.next_in = reader->chunk_in;
    reader->stream.avail_in += read;
}

/// consumes as much of avail_in as possible, writing it to chunk_out
/// must only be called once all of chunk_out has been fully read (read_out must be == size-stream.avail_out).
static void decompress_into_chunk(CompressedReader *reader) {
    assert(
        reader->read_out == reader->chunk_size - reader->stream.avail_out
    );
    reader->stream.avail_out = reader->chunk_size;
    reader->stream.next_out = reader->chunk_out;
    reader->read_out = 0;
    const int ret = inflate(&reader->stream, Z_SYNC_FLUSH);
    switch (ret) {
        case Z_NEED_DICT:
        case Z_DATA_ERROR:
        case Z_STREAM_ERROR:
            inflateEnd(&reader->stream);
            reader->is_err = true;
            reader->error_msg = "Couldn't decompress. Incomplete or corrupt zlib compressed data.";
            return;
        case Z_MEM_ERROR:
            inflateEnd(&reader->stream);
            reader->is_err = true;
            reader->error_msg = "Couldn't allocate memory.";
            return;
        default:;
    }
}

/// will read n compressed bytes, and error if it is unable to do so.
static char *read_bytes(CompressedReader *reader, const size_t n) {
    char *res = malloc(n+1);
    res[n] = 0;
    for (size_t i = 0; i < n; i++) {
        if (reader->read_out == reader->chunk_size - reader->stream.avail_out) {
            // we ran out of chunk_out
            if (reader->stream_over)
            {
                reader->is_err = true;
                reader->error_msg = "Stream ended prematurely.";
                return res;
            }
            read_file_into_chunk(reader);
            decompress_into_chunk(reader);
        }

        res[i] = (char)reader->chunk_out[reader->read_out];
        reader->read_out++;
    }
    if (reader->stream_over)
    {
        inflateEnd(&reader->stream);
    }
    reader->crc = crc32(reader->crc, (Bytef *)res, n);
    return res;
}


// names based on the python "struct" codes
static CompressedReader *new_compressed_reader(const char *path) {
    CompressedReader *res = (CompressedReader *)calloc(sizeof(CompressedReader), 1);
    res->is_err = false;
    res->stream_over = false;
    FILE *file = fopen(path, "rb");
    if (file == NULL) {
        res->is_err = true;
        res->error_msg = "File not found.";
        return res;
    }
    res->file = file;
    res->crc = crc32(0L, Z_NULL, 0);

    const size_t uncompressed_header_size = 8;
    res->header = calloc(uncompressed_header_size, sizeof(char));
    if (fread(res->header, 1, uncompressed_header_size, file) < uncompressed_header_size) {
        res->is_err = true;
        res->error_msg = "File too short, can't read header.";
        return res;
    }
    res->crc = crc32(res->crc, (Bytef *)res->header, uncompressed_header_size);

    res->chunk_size = 1 << 20; // 1 MB
    res->chunk_in = (Bytef *)malloc(res->chunk_size);
    res->chunk_out = (Bytef *)malloc(res->chunk_size);

    res->stream.zalloc = Z_NULL;
    res->stream.zfree = Z_NULL;
    res->stream.opaque = Z_NULL;
    res->stream.avail_in = 0;
    res->stream.next_in = NULL;
    res->stream.avail_out = res->chunk_size;
    res->stream.next_out = res->chunk_out;

    read_file_into_chunk(res);

    if (inflateInit(&res->stream) != Z_OK) {
        res->is_err = true;
        res->error_msg = "Couldn't initialize zlib stream.";
        return res;
    }

    decompress_into_chunk(res);

    return res;
}

static void free_compressed_reader(CompressedReader *reader) {
    if (reader->chunk_out != NULL)
        free(reader->chunk_out);
    if (reader->chunk_in != NULL)
        free(reader->chunk_in);
    if (reader->file != NULL) {
        fclose(reader->file);
        reader->file = NULL;
    }
}

static char *read_s(CompressedReader *reader) {
    char *n_bytes_buf = read_bytes(reader, 1);
    const unsigned char n_bytes = *n_bytes_buf;
    free(n_bytes_buf);

    char *res = read_bytes(reader, n_bytes);
    return res;
}
static size_t read_Q(CompressedReader *reader) {
    char *n_bytes_buf = read_bytes(reader, 8);
    // we assume little endian machine
    const size_t res = *(size_t *)n_bytes_buf;
    free(n_bytes_buf);
    return res;
}
static uint32_t read_I(CompressedReader *reader) {
    char *n_bytes_buf = read_bytes(reader, 4);
    // we assume little endian machine
    const uint32_t res = *(uint32_t *)n_bytes_buf;
    free(n_bytes_buf);
    return res;
}
static float read_f(CompressedReader *reader) {
    char *n_bytes_buf = read_bytes(reader, 4);
    // we assume little endian machine
    const float res = *(float *)n_bytes_buf;
    free(n_bytes_buf);
    return res;
}
static double read_d(CompressedReader *reader) {
    char *n_bytes_buf = read_bytes(reader, 8);
    // we assume little endian machine
    const double res = *(double *)n_bytes_buf;
    free(n_bytes_buf);
    return res;
}
static char *read_Is(CompressedReader *reader, const size_t n) {
    size_t cap = n * 8;
    size_t len = 0;
    char *res = malloc(cap+1);
    for (size_t i = 0; i < n; i++) {
        const uint32_t I = read_I(reader);
        int extra = snprintf(&res[len], cap-len, "%u ", I);
        assert(extra >= 0);
        while (len + extra > cap)
        {
            cap *= 2;
            res = realloc(res, cap+1);
            assert(res != NULL);
        }
        extra = snprintf(&res[len], cap-len, "%u ", I);
        assert(extra >= 0);
        assert(len + extra <= cap);
        len += extra;
    }
    res = realloc(res, len+1);
    assert(res != NULL);
    res[len] = 0;
    return res;
}
static char *read_fs(CompressedReader *reader, const size_t n) {
    size_t cap = n * 8;
    size_t len = 0;
    char *res = malloc(cap+1);
    for (size_t i = 0; i < n; i++) {
        const float f = read_f(reader);
        int extra = snprintf(&res[len], cap-len, "%.1f ", f);
        assert(extra >= 0);
        while (len + extra > cap)
        {
            cap *= 2;
            res = realloc(res, cap+1);
            assert(res != NULL);
        }
        extra = snprintf(&res[len], cap-len, "%.1f ", f);
        assert(extra >= 0);
        assert(len + extra <= cap);
        len += extra;
    }
    res = realloc(res, len+1);
    assert(res != NULL);
    res[len] = 0;
    return res;
}
static char *read_ss(CompressedReader *reader, const size_t n) {
    size_t cap = n * 8;
    size_t len = 0;
    char *res = malloc(cap+1);
    for (size_t i = 0; i < n; i++) {
        char *s = read_s(reader);
        const size_t extra = strlen(s);
        while (len + extra + 1 > cap)
        {
            cap *= 2;
            res = realloc(res, cap+1);
            assert(res != NULL);
        }
        memcpy(&res[len], s, extra);
        res[len+extra] = ' ';
        len += extra + 1;
        free(s);
    }
    res = realloc(res, len+1);
    assert(res != NULL);
    res[len] = 0;
    return res;
}

typedef struct bond {
    uint32_t bi;
    uint32_t bj;
} Bond;

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

static char *read_bonds(CompressedReader *reader, size_t n_atoms) {
    const size_t n_bonds = read_Q(reader);

    // sorted list of bonds
    size_t sorted_len = 0;
    Bond *sorted_bonds = calloc(sizeof(struct bond), 2 * n_bonds);
    
    for (size_t i = 0; i < n_bonds; i++)
    {
        Bond b = { read_I(reader), read_I(reader) };
        if (b.bi == b.bj) {
            // self bonding
            // forbidden by the file format, but we silently ignore it
            continue;
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

static bool read_crc(CompressedReader *reader) {
    const uint32_t reference = reader->crc;
    const uint32_t crc = read_I(reader);
    reader->crc = crc32(0L, Z_NULL, 0);;
    return crc == reference;
}

/* Global state */
// sadly required, as the traces have to be reading information from somewhere

typedef struct {
    uint32_t frame_number;
    uint32_t n_atoms;
    size_t sim_step;
    double sim_time_ns;
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
    bool is_err;
    char *error_msg;

    // metadata of how it was made
    char *path;
    int molid;
    int remove_pbc_crossing;
    // data loaded into memory
    // per simulation data

    // per frame data
    // number of frames loaded
    size_t n_frames;
    // current capacity of frames
    size_t cap_frames;
    TopTrajFrame **frames;
} TopTrajData;

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

static TopTrajData *load_toptraj(const char *path, int molid, int pbc) {
    /// Reads .toptraj file at path and loads it into a dynamically allocated
    /// object, which it returns.
    /// While doing so, prints a progress bar on STDOUT.
    /// Prints error message to STDERR and returns NULL if there is any error.
    TopTrajData *res = (TopTrajData *)calloc(sizeof(TopTrajData), 1);
    res->is_err = false;

    res->molid = molid;
    res->remove_pbc_crossing = pbc;
    // clone path, as we do not own it
    res->path = malloc(strlen(path) + 1);

    strcpy(res->path, path);

#define CHECK_ERR \
if (reader->is_err) { \
    res->is_err = true; \
    res->error_msg = reader->error_msg; \
    free_compressed_reader(reader); \
    return res; \
}
    CompressedReader *reader = new_compressed_reader(path); CHECK_ERR
    // 1.0 magic number
    if (strncmp(reader->header, "\xc0TOPTR\x01\0", 8) == 0) {
        printf("TopTraj file version is 1.0\n");

        char *name = read_s(reader); CHECK_ERR
        printf("Simulation name: %s\n", name);
        const uint32_t n_initial_molecules = read_I(reader); CHECK_ERR
        printf("Number of initial molecule types: %u\n", n_initial_molecules);

        for (uint32_t i = 0; i < n_initial_molecules; i++)
        {
            char *molname = read_s(reader);
            const uint32_t molcount = read_I(reader);
            printf("Mol name: %s\n", molname);
            printf("Mol count: %u\n", molcount);
        }

        if (!read_crc(reader)) {
            free_compressed_reader(reader);
            res->is_err = true;
            res->error_msg = "CRC mismatch in header! .toptraj file is likely corrupt.";
            return res;
        }

        res->cap_frames = 1000;
        res->frames = calloc(sizeof(TopTrajFrame), res->cap_frames);
        res->n_frames = 0;

        size_t frame_index = 0;
        printf("\n");
        while (!reader->stream_over) {
            TopTrajFrame *frame = calloc(sizeof(TopTrajFrame), 1);
            frame->frame_number = read_I(reader);
            if (frame->frame_number != frame_index) {
                // probably junk at the end of the file from a sim that
                // was cancelled
                printf("\nFrame %lu wrong frame index, stopping.", frame_index + 1);
                free(frame);
                break;
            }
            frame->n_atoms = read_I(reader);
            frame->sim_step = read_Q(reader);
            frame->sim_time_ns = read_d(reader);

            frame->names = read_ss(reader, frame->n_atoms);
            frame->resnames = read_ss(reader, frame->n_atoms);
            frame->resids = read_Is(reader, frame->n_atoms);
            frame->types = read_ss(reader, frame->n_atoms);
            frame->charges = read_fs(reader, frame->n_atoms);
            frame->masses = read_fs(reader, frame->n_atoms);
            frame->bonds = read_bonds(reader, frame->n_atoms);
            if (!read_crc(reader)) {
                free_frame(frame);
                printf("Frame %lu CRC mismatch, stopping.", frame_index + 1);
                break;
            }
            printf("\rFrame read successfully: %lu", frame_index + 1);
            frame_index++;
            res->frames[res->n_frames] = frame;
            res->n_frames++;
            if (res->n_frames >= res->cap_frames)
            {
                res->cap_frames *= 2;
                res->frames = realloc(res->frames, sizeof(TopTrajFrame) * res->cap_frames);
                assert(res->frames != NULL);
            }

        }
        printf("\nRead %lu frames.\n", frame_index);

        free_compressed_reader(reader);
        return res;

    } else {
        free_compressed_reader(reader);
        res->is_err = true;
        res->error_msg = "File is not a valid toptraj file, or the version is unknown.";
        return res;
    }
}



static char *on_frame_change(
    ClientData data, // TopTrajData *
    Tcl_Interp *interp,
    const char *name1, // "vmd_frame"
    const char *name2,
    int flags
) {
    TopTrajData *toptraj = (TopTrajData *)data;

    // GET CURRENT FRAME
    char molid[64];
    snprintf(molid, 64, "%i", toptraj->molid);

    int res = Tcl_VarEval(interp, "molinfo ", molid, " get frame", NULL);
    if (res != TCL_OK) {
        printf("TOPTRAJ FATAL: Trace failed at molinfo get frame: %s\n", Tcl_GetStringResult(interp));
        return NULL;
    }
    int64_t c_frame = atol(Tcl_GetStringResult(interp));
    if (c_frame < 0 || c_frame > toptraj->n_frames) {
        printf("TOPTRAJ FATAL: Frame %li is out of range for the loaded .toptraj.\n", c_frame);
        return NULL;
    }

    char frame[64];
    snprintf(frame, 64, "%li", c_frame);

    // ATOMSELECT
    res = Tcl_VarEval(interp, "atomselect ", molid, " all frame ", frame, NULL);
    if (res != TCL_OK) {
        printf("TOPTRAJ FATAL: Trace failed at atomselect: %s\n", Tcl_GetStringResult(interp));
        return NULL;
    }
    const char *sel = Tcl_GetStringResult(interp);

    // N_ATOMS
    res = Tcl_VarEval(interp, sel, " num", NULL);
    if (res != TCL_OK) {
        printf("TOPTRAJ FATAL: Trace failed at n_atoms: %s\n", Tcl_GetStringResult(interp));
        return NULL;
    }
    int64_t n_atoms = atol(Tcl_GetStringResult(interp));
    if (n_atoms != toptraj->frames[c_frame]->n_atoms) {
        printf("TOPTRAJ FATAL: number of atoms does not match .toptraj.\n");
        return NULL;
    }

    // SETBONDS
    res = Tcl_VarEval(interp, sel, " setbonds ", toptraj->frames[c_frame]->bonds, NULL);
    if (res != TCL_OK) {
        printf("TOPTRAJ FATAL: Trace failed at setbonds: %s\n", Tcl_GetStringResult(interp));
        return NULL;
    }

    res = Tcl_VarEval(interp, sel, " delete", NULL);
    if (res != TCL_OK) {
        printf("TOPTRAJ FATAL: Trace failed at delete: %s\n", Tcl_GetStringResult(interp));
        return NULL;
    }

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
    TopTrajData *data = load_toptraj(argv[1], molid, pbc);
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

static int unload_toptraj_cmd(
    ClientData _data, // NULL
    Tcl_Interp *interp,
    int argc, char const *argv[]
) {
    (void)_data; // supress unused warning
    Tcl_SetResult(interp, "Unimplemented", TCL_STATIC); // TODO
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
