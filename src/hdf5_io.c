#include <hdf5.h>
#include <quadmath.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

static int write_string_dataset(hid_t fid, const char *name, const char *value) {
  hid_t dtype = -1, space = -1, dset = -1;
  herr_t st;
  size_t len = strlen(value);

  dtype = H5Tcopy(H5T_C_S1);
  if (dtype < 0) goto fail;
  if (H5Tset_size(dtype, len) < 0) goto fail;
  if (H5Tset_strpad(dtype, H5T_STR_SPACEPAD) < 0) goto fail;

  space = H5Screate(H5S_SCALAR);
  if (space < 0) goto fail;

  dset = H5Dcreate2(fid, name, dtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dset < 0) goto fail;

  st = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
  if (st < 0) goto fail;

  H5Dclose(dset);
  H5Sclose(space);
  H5Tclose(dtype);
  return 1;

fail:
  if (dset >= 0) H5Dclose(dset);
  if (space >= 0) H5Sclose(space);
  if (dtype >= 0) H5Tclose(dtype);
  return 0;
}

int hdf5_write_string_pair(const char *file,
                           const char *name1, const char *value1,
                           const char *name2, const char *value2) {
  hid_t fid = -1;
  int ok = 0;

  fid = H5Fcreate(file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (fid < 0) return 0;

  ok = write_string_dataset(fid, name1, value1);
  if (ok) ok = write_string_dataset(fid, name2, value2);

  H5Fclose(fid);
  return ok;
}

static void format_float128_strings(char *out, const __float128 *vals, size_t n) {
  size_t i;
  char tmp[64];

  for (i = 0; i < n; ++i) {
    memset(out + 64*i, ' ', 64);
    memset(tmp, 0, sizeof(tmp));
    quadmath_snprintf(tmp, sizeof(tmp), "%+-#46.36QE", vals[i]);
    memcpy(out + 64*i, tmp, strlen(tmp));
  }
}

static int write_float128_string_array(hid_t fid, const char *name,
                                       const __float128 *vals, int rank,
                                       const hsize_t *dims) {
  hid_t dtype = -1, space = -1, dset = -1;
  hsize_t count = 1;
  char *buf = NULL;
  int i;
  herr_t st;

  for (i = 0; i < rank; ++i) count *= dims[i];

  buf = (char *)malloc((size_t)count * 64);
  if (buf == NULL) goto fail;
  format_float128_strings(buf, vals, (size_t)count);

  dtype = H5Tcopy(H5T_C_S1);
  if (dtype < 0) goto fail;
  if (H5Tset_size(dtype, 64) < 0) goto fail;
  if (H5Tset_strpad(dtype, H5T_STR_SPACEPAD) < 0) goto fail;

  space = H5Screate_simple(rank, dims, NULL);
  if (space < 0) goto fail;

  dset = H5Dcreate2(fid, name, dtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dset < 0) goto fail;

  st = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  if (st < 0) goto fail;

  free(buf);
  H5Dclose(dset);
  H5Sclose(space);
  H5Tclose(dtype);
  return 1;

fail:
  if (buf != NULL) free(buf);
  if (dset >= 0) H5Dclose(dset);
  if (space >= 0) H5Sclose(space);
  if (dtype >= 0) H5Tclose(dtype);
  return 0;
}

int hdf5_write_legeexps_r128(const char *file, int64_t n,
                             const __float128 *x,
                             const __float128 *u,
                             const __float128 *v,
                             const __float128 *whts) {
  hid_t fid = -1;
  hsize_t dim1[1], dim2[2];
  int ok = 0;

  fid = H5Fcreate(file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (fid < 0) return 0;

  dim1[0] = (hsize_t)n;
  dim2[0] = (hsize_t)n;
  dim2[1] = (hsize_t)n;

  ok = write_float128_string_array(fid, "/x", x, 1, dim1);
  if (ok) ok = write_float128_string_array(fid, "/u", u, 2, dim2);
  if (ok) ok = write_float128_string_array(fid, "/v", v, 2, dim2);
  if (ok) ok = write_float128_string_array(fid, "/whts", whts, 1, dim1);

  H5Fclose(fid);
  return ok;
}

int hdf5_write_real128_matrix(const char *file, const char *name,
                              int64_t n1, int64_t n2,
                              const __float128 *vals) {
  hid_t fid = -1;
  hsize_t dims[2];
  int ok;

  fid = H5Fcreate(file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (fid < 0) return 0;

  dims[0] = (hsize_t)n1;
  dims[1] = (hsize_t)n2;
  ok = write_float128_string_array(fid, name, vals, 2, dims);

  H5Fclose(fid);
  return ok;
}

int hdf5_write_two_real128_matrices(const char *file,
                                    const char *name1, const __float128 *vals1,
                                    const char *name2, const __float128 *vals2,
                                    int64_t n1, int64_t n2) {
  hid_t fid = -1;
  hsize_t dims[2];
  int ok;

  fid = H5Fcreate(file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (fid < 0) return 0;

  dims[0] = (hsize_t)n1;
  dims[1] = (hsize_t)n2;
  ok = write_float128_string_array(fid, name1, vals1, 2, dims);
  if (ok) ok = write_float128_string_array(fid, name2, vals2, 2, dims);

  H5Fclose(fid);
  return ok;
}

/* ------------------------------------------------------------------
 * Generic write-one-array helper. If append_flag == 0, the file is
 * truncated (created fresh). If append_flag != 0, the file is opened
 * R/W and the dataset is added (caller must ensure unique name).
 *
 * `dims` is a length-`rank` array of 64-bit dataset shape (column-
 * major from Fortran callers; same layout as the matrix writers).
 * ------------------------------------------------------------------ */
int hdf5_write_real128_array(const char *file, const char *name,
                             int rank, const int64_t *dims,
                             const __float128 *vals, int append_flag) {
  hid_t fid = -1;
  hsize_t hdims[8];
  int i, ok;

  if (rank < 1 || rank > 8) return 0;
  for (i = 0; i < rank; ++i) hdims[i] = (hsize_t)dims[i];

  if (append_flag == 0) {
    fid = H5Fcreate(file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  } else {
    fid = H5Fopen(file, H5F_ACC_RDWR, H5P_DEFAULT);
  }
  if (fid < 0) return 0;

  ok = write_float128_string_array(fid, name, vals, rank, hdims);

  H5Fclose(fid);
  return ok;
}

/* ------------------------------------------------------------------
 * Read a r128 dataset back into __float128 array.
 *
 * Storage layout matches the writers (length-64 fixed-string per
 * value, formatted via quadmath %+-#46.36QE). Caller passes
 * expected `rank` and `dims`; if the dataset's actual rank/dims
 * differ, the function fails (returns 0).
 *
 * `vals` must be pre-allocated by the caller to fit the full array
 * (product of dims) of __float128.
 * ------------------------------------------------------------------ */
/* ------------------------------------------------------------------
 * Query the rank/dims of a r128 dataset without reading data.
 *
 * Caller passes the expected `rank`. On success, fills `dims` with the
 * dataset's actual extent and returns 1; on rank mismatch or any HDF5
 * failure, returns 0. Used by the per-triangle r128 mex routines so
 * they can size their full-mesh buffers from the file metadata instead
 * of taking an `ntri_full` arg from MATLAB.
 * ------------------------------------------------------------------ */
int hdf5_get_dims_r128(const char *file, const char *name,
                       int rank, int64_t *dims) {
  hid_t fid = -1, dset = -1, space = -1;
  hsize_t actual_dims[8];
  int actual_rank, i, ret = 0;

  if (rank < 1 || rank > 8) return 0;

  fid = H5Fopen(file, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (fid < 0) goto fail;

  dset = H5Dopen2(fid, name, H5P_DEFAULT);
  if (dset < 0) goto fail;

  space = H5Dget_space(dset);
  if (space < 0) goto fail;

  actual_rank = H5Sget_simple_extent_ndims(space);
  if (actual_rank != rank) goto fail;
  if (H5Sget_simple_extent_dims(space, actual_dims, NULL) < 0) goto fail;

  for (i = 0; i < rank; ++i) dims[i] = (int64_t)actual_dims[i];
  ret = 1;

fail:
  if (space >= 0) H5Sclose(space);
  if (dset >= 0) H5Dclose(dset);
  if (fid >= 0) H5Fclose(fid);
  return ret;
}

int hdf5_read_real128_array(const char *file, const char *name,
                            int rank, const int64_t *dims,
                            __float128 *vals) {
  hid_t fid = -1, dset = -1, space = -1, dtype = -1;
  hsize_t actual_dims[8];
  hsize_t count = 1;
  int actual_rank;
  char *buf = NULL;
  size_t i, str_size;
  int ret = 0;

  if (rank < 1 || rank > 8) return 0;
  for (i = 0; (int)i < rank; ++i) count *= (hsize_t)dims[i];

  fid = H5Fopen(file, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (fid < 0) goto fail;

  dset = H5Dopen2(fid, name, H5P_DEFAULT);
  if (dset < 0) goto fail;

  space = H5Dget_space(dset);
  if (space < 0) goto fail;

  actual_rank = H5Sget_simple_extent_ndims(space);
  if (actual_rank != rank) goto fail;
  if (H5Sget_simple_extent_dims(space, actual_dims, NULL) < 0) goto fail;
  for (i = 0; (int)i < rank; ++i) {
    if (actual_dims[i] != (hsize_t)dims[i]) goto fail;
  }

  dtype = H5Dget_type(dset);
  if (dtype < 0) goto fail;
  str_size = H5Tget_size(dtype);
  if (str_size != 64) goto fail;

  buf = (char *)malloc((size_t)count * 64);
  if (buf == NULL) goto fail;

  if (H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf) < 0) goto fail;

  /* Parse each 64-byte slot back to __float128. */
  for (i = 0; i < (size_t)count; ++i) {
    char tmp[65];
    memcpy(tmp, buf + 64*i, 64);
    tmp[64] = '\0';
    vals[i] = strtoflt128(tmp, NULL);
  }

  ret = 1;

fail:
  if (buf != NULL) free(buf);
  if (dtype >= 0) H5Tclose(dtype);
  if (space >= 0) H5Sclose(space);
  if (dset >= 0) H5Dclose(dset);
  if (fid >= 0) H5Fclose(fid);
  return ret;
}
