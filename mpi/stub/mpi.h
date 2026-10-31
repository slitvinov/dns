/* mpi.h: single-process stand-in for the MPI calls of boxsymmp (np 1 only; build with -Istub where no MPI exists) */
#pragma once
#include <string.h>
#include <stdlib.h>
#include <time.h>
typedef int MPI_Comm;
typedef int MPI_Datatype;
typedef int MPI_Op;
typedef int MPI_Info;
enum { MPI_COMM_WORLD = 0, MPI_COMM_TYPE_SHARED = 0, MPI_INFO_NULL = 0, MPI_SUM = 0, MPI_MAX = 1 };
enum { MPI_BYTE = 1, MPI_INT = 4, MPI_FLOAT = 5, MPI_DOUBLE = 8 }; /* value = size in bytes (MPI_FLOAT: 4, see sz) */
#define MPI_IN_PLACE ((void *)1)
static int mpistub_tsize[64];
static inline size_t mpistub_sz(MPI_Datatype t) { return t == MPI_FLOAT ? 4 : t < 64 && mpistub_tsize[t] ? mpistub_tsize[t] : (size_t)t; }
static inline int MPI_Init(int *, char ***) { return 0; }
static inline int MPI_Finalize(void) { return 0; }
static inline int MPI_Comm_rank(MPI_Comm, int *r) { *r = 0; return 0; }
static inline int MPI_Comm_size(MPI_Comm, int *s) { *s = 1; return 0; }
static inline int MPI_Comm_split_type(MPI_Comm, int, int, MPI_Info, MPI_Comm *c) { *c = 0; return 0; }
static inline int MPI_Abort(MPI_Comm, int e) { exit(e); }
static inline int MPI_Barrier(MPI_Comm) { return 0; }
static inline int MPI_Bcast(void *, int, MPI_Datatype, int, MPI_Comm) { return 0; }
static int mpistub_next = 16;
static inline int MPI_Type_contiguous(int n, MPI_Datatype t, MPI_Datatype *o) { *o = mpistub_next++; mpistub_tsize[*o] = n * mpistub_sz(t); return 0; }
static inline int MPI_Type_commit(MPI_Datatype *) { return 0; }
static inline double MPI_Wtime(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
static inline int MPI_Allreduce(const void *s, void *r, int n, MPI_Datatype t, MPI_Op, MPI_Comm) {
  if (s != MPI_IN_PLACE) memcpy(r, s, n * mpistub_sz(t)); return 0; }
static inline int MPI_Gatherv(const void *s, int n, MPI_Datatype t, void *r, const int *, const int *d, MPI_Datatype, int, MPI_Comm) {
  memcpy((char *)r + (size_t)d[0] * mpistub_sz(t), s, n * mpistub_sz(t)); return 0; }
static inline int MPI_Alltoallv(const void *, const int *, const int *, MPI_Datatype, void *, const int *, const int *, MPI_Datatype, MPI_Comm) {
  abort(); } /* BOXSYM_MPI=1 is not available with the stub */
