#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <tcl/tcl.h>
// switch which lines are commented on e.g. fedora
//#include <tcl.h>
#include <tcl/tclDecls.h>
//#include <tclDecls.h>
#include <zconf.h>
#include <zlib.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

ssize_t load_file(const char *path, char **out) {
  // Loads a file at path, returns its length and sets the (char *) buffer out
  // to point to the contents. If there is an error, it returns -1 and the out
  // parameter will point to the error message (in static memory)
  #define ERROR(msg) \
    *out = msg; \
    return -1;

  FILE *file = fopen(path, "rb");
  if (!file) {
    ERROR("File not found.");
  }
  // get file length by going to the end
  fseek(file, 0L, SEEK_END);
  ssize_t file_len = ftell(file);
  if (file_len > INT32_MAX) {
    ERROR("Compressed file can be at most 2GB large.");
  }
  // seek back to character 0
  fseek(file, 0L, SEEK_SET);

  // alloc memory for the file and read the file
  *out = malloc(file_len);
  if (!*out) {
    ERROR("Memory allocation failed.");
  }
  fread(*out, 1, file_len, file);
  if (ferror(file)) {
    free(*out);
    *out = NULL;
    ERROR("Error reading file.");
  }
  fclose(file);

  return file_len;
  #undef ERROR
}

ssize_t decompress(ssize_t file_len, Bytef *compressed, char **out) {
  // Returns the number of bytes decompressed, and sets out to point to the
  // output buffer. In should be cleaned up by caller.
  // If there is an error, returns -1 and sets out to point to the error
  // message (in static memory)
  // decompress
  // buffer - 128 MB
  #define ERROR(msg) \
    *out = msg; \
    return -1;
  #define SIZE 134217728
  Bytef *chunk = (Bytef *)malloc(SIZE);
  // output
  *out = malloc(SIZE);
  ssize_t out_cap = SIZE;
  ssize_t out_len = 0;
  // zlib stream
  z_stream strm;
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  strm.avail_in = file_len;
  strm.next_in = compressed;

  if (inflateInit(&strm) != Z_OK) {
    ERROR("Couldn't initialize zlib stream.");
  }

  int ret;
  do {
    // decompress a chunk's worth
    strm.avail_out = SIZE;
    strm.next_out = &chunk[0];
    ret = inflate(&strm, Z_NO_FLUSH);
    switch (ret) {
      case Z_NEED_DICT:
      case Z_DATA_ERROR:
      case Z_STREAM_ERROR:
        inflateEnd(&strm);
        free(chunk);
        free(*out);
        ERROR("Invalid, incomplete or corrupt zlib compressed data.")
      case Z_MEM_ERROR:
        inflateEnd(&strm);
        free(chunk);
        free(*out);
        ERROR("Error decompressing, out of memory.")
      default:;
    }
    // copy to output buffer
    ssize_t new = SIZE - strm.avail_out;
    if (out_len + new > out_cap) {
      out_cap += SIZE;
      if (!(*out = realloc(*out, out_cap))) {
        inflateEnd(&strm);
        free(chunk);
        ERROR("Error allocating memory.");
      }
    }
    memcpy(&(*out[out_len]), &chunk[0], new);
    out_len += new;
  } while (ret != Z_STREAM_END);
  // cleanup
  inflateEnd(&strm);
  free(chunk);
  if (!(*out = realloc(*out, out_len))) {
    ERROR("Error shrinking allocation.");
  }
  return out_len;
  #undef ERROR
}

typedef struct bonds_info {
  bool success;
  char *msg;
  ssize_t n_atoms;
  ssize_t n_frames;
  char **frames;
  ssize_t *n_chars;
  char *body;
} BondsInfo;

BondsInfo parse(ssize_t n_atoms, ssize_t n_frames, char *source) {
  BondsInfo res;
  res.success = true;
  res.msg = NULL;
  res.n_atoms = n_atoms;
  res.n_frames = n_frames;
  res.body = source;
  #define ERROR(message) \
    res.success = false; \
    res.msg = message;
  res.frames = malloc(sizeof(char *) * n_frames);
  res.n_chars = malloc(sizeof(ssize_t) * n_frames);
  if (!res.frames) {
    ERROR("Memory allocation failed.");
  }
  char *c = source;
  ssize_t atom_index = 0;
  ssize_t frame_index = -1;
  while (*c != '\0') {
    if (frame_index >= n_frames) {
      break;
    }
    switch (*c) {
      case '{':
        if (atom_index == 0) {
          // start of a frame
          frame_index++;
          res.frames[frame_index] = c;
        }
        break;
      case '}':
        atom_index++;
        if (atom_index == n_atoms) {
          // end of a frame
          res.n_chars[frame_index] = c - res.frames[frame_index] + 1;
          atom_index = 0;
        }
        break;
      default:;
    }
    c++;
  }
  if (atom_index != 0) {
    ERROR("Provided n_atoms is likely wrong (not a divisor of the number of bonds in file).")
  }
  if (frame_index + 1 != n_frames) {
    printf("Frames parsed: %li, frames expected: %li. Is this the right trajectory file?\n", frame_index+1, n_frames);
    ERROR("Wrong frames length.")
  }
  return res;
  #undef ERROR
}

BondsInfo global_bonds;

void free_buf(char *buf) {
  free(buf);
}

static int Bonds_Frame_Cmd(ClientData cdata, Tcl_Interp *interp, int argc, char const *argv[]) {
  #define ERROR(msg) \
  Tcl_SetResult(interp, msg, TCL_STATIC); \
  return TCL_ERROR;
  if (argc < 2) {
    ERROR("Usage: bonds_frame [frame index].")
  }
  if (*argv[1] < '0' | *argv[1] > '9') {
    ERROR("Please provide a number argument to bonds_frame.")
  }
  ssize_t frame_index = atol(argv[1]);
  if (!global_bonds.body || !global_bonds.frames || !global_bonds.n_chars) {
    ERROR("No bonds loaded. Please load them with bonds_load first.")
  }
  if (frame_index < 0 || frame_index >= global_bonds.n_frames) {
    ERROR("Index for bonds_frame out of range.")
  }
  char *buf = calloc(global_bonds.n_chars[frame_index] + 1, sizeof(char));
  memcpy(buf, global_bonds.frames[frame_index], global_bonds.n_chars[frame_index]);
  Tcl_SetResult(interp, buf, free_buf);
  return TCL_OK;
  #undef ERROR
}


void free_globals() {
  if (global_bonds.body != NULL) {
    free(global_bonds.body);
    global_bonds.body = NULL;
  }
  if (global_bonds.frames != NULL) {
    free(global_bonds.frames);
    global_bonds.frames = NULL;
  }
  if (global_bonds.n_chars != NULL) {
    free(global_bonds.n_chars);
    global_bonds.n_chars = NULL;
  }
}

static int Bonds_Free_Cmd(ClientData cdata, Tcl_Interp *interp, int argc, char const *argv[]) {
  free_globals();
  return TCL_OK;
}

static int Load_Bonds_Cmd(ClientData cdata, Tcl_Interp *interp, int argc, char const *argv[]) {
  // the "bonds_load" command, as exposed to Tcl
  // it loads the bond info from <path> to a global variable in this file
  // later, bonds_frame can be used to get a single frame
  // bonds_free can be called to free the global bonds info
  #define ERROR(msg) \
    Tcl_SetResult(interp, msg, TCL_STATIC); \
    return TCL_ERROR;

  if (argc < 4) {
    ERROR("Usage: bonds_load [filename] [n atoms] [n frames].");
  }
  const char *path = argv[1];
  ssize_t n_atoms = atol(argv[2]);
  ssize_t n_frames = atol(argv[3]);
  if (n_atoms <= 0 || n_frames <= 0) {
    ERROR("n_atoms and n_frames (args 2 and 3) have to be a positive non-zero number.");
  }

  // open file
  Bytef *compressed;
  assert (sizeof(char) == sizeof(Bytef));
  ssize_t file_len = load_file(path, (char **)&compressed);
  if (file_len == -1) {
    ERROR((char *)compressed);
  }
  printf("File %s successfully read into memory (%li bytes on disk).\n", path, file_len);

  // decompress file, TODO stream the input file to reduce memory usage?
  char *decompressed;
  ssize_t buf_len;
  buf_len = decompress(file_len, compressed, &decompressed);
  free(compressed);
  if (buf_len == -1) {
    ERROR((char *)decompressed);
  }
  printf("File %s successfully decompressed (%li bytes).\n", path, buf_len);

  // parse file
  BondsInfo new_bonds_info = parse(n_atoms, n_frames, decompressed);
  if (!new_bonds_info.success) {
    free(decompressed);
    ERROR(new_bonds_info.msg);
  }
  free_globals();
  global_bonds = new_bonds_info;
  printf("File %s successfully parsed (%li frames).\n", path, global_bonds.n_frames);
  return TCL_OK;
  #undef ERROR
}

int DLLEXPORT Daemon_Init(Tcl_Interp *interp) {
  if (Tcl_InitStubs(interp, TCL_VERSION, 0) == NULL) {
    return TCL_ERROR;
  }
  // TODO tcl_pkgprovide ?
  Tcl_CreateCommand(interp, "bonds_load", Load_Bonds_Cmd, NULL, NULL);
  Tcl_CreateCommand(interp, "bonds_frame", Bonds_Frame_Cmd, NULL, NULL);
  Tcl_CreateCommand(interp, "bonds_free", Bonds_Free_Cmd, NULL, NULL);
  memset(&global_bonds, 0, sizeof(BondsInfo));
  return TCL_OK;
}
