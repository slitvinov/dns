/* nccl.h: single-process stand-in for the NCCL calls of boxsymmp (np 1: no peers, so send/recv never run) */
#pragma once
#include <cuda_runtime.h>
typedef int ncclComm_t;
typedef enum { ncclSuccess = 0 } ncclResult_t;
typedef enum { ncclFloat = 7, ncclDouble = 8 } ncclDataType_t;
typedef struct { char internal[128]; } ncclUniqueId;
static inline ncclResult_t ncclGetUniqueId(ncclUniqueId *) { return ncclSuccess; }
static inline ncclResult_t ncclCommInitRank(ncclComm_t *c, int, ncclUniqueId, int) { *c = 0; return ncclSuccess; }
static inline ncclResult_t ncclCommDestroy(ncclComm_t) { return ncclSuccess; }
static inline const char *ncclGetErrorString(ncclResult_t) { return "nccl stub"; }
static inline ncclResult_t ncclGroupStart(void) { return ncclSuccess; }
static inline ncclResult_t ncclGroupEnd(void) { return ncclSuccess; }
static inline ncclResult_t ncclSend(const void *, size_t, ncclDataType_t, int, ncclComm_t, cudaStream_t) { return ncclSuccess; }
static inline ncclResult_t ncclRecv(void *, size_t, ncclDataType_t, int, ncclComm_t, cudaStream_t) { return ncclSuccess; }
