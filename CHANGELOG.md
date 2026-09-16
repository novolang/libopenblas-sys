# Changelog

All notable changes to libopenblas-sys are recorded here. The format
is [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
package follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
with the pre-1.0 rule that a breaking change bumps the MINOR number.

## 0.1.0 — 2026-09-16

The first release: twenty-seven entry points of libopenblas through the
CBLAS interface, one `@ffi` declaration each, and no logic.

### Added

- `libopenblas` — the whole surface, in four groups.
  - Level 1, vector against vector: `cblas_ddot`, `cblas_dnrm2`,
    `cblas_dasum`, `cblas_idamax`, `cblas_daxpy`, `cblas_dscal`,
    `cblas_dcopy`, `cblas_dswap`, `cblas_drot` and `cblas_drotg`.
  - Level 2, matrix against vector: `cblas_dgemv`, `cblas_dger`,
    `cblas_dsymv`, `cblas_dtrmv` and `cblas_dtrsv`.
  - Level 3, matrix against matrix: `cblas_dgemm`, `cblas_dsymm`,
    `cblas_dsyrk`, `cblas_dsyr2k`, `cblas_dtrmm` and `cblas_dtrsm`.
  - The OpenBLAS runtime: `openblas_set_num_threads`,
    `openblas_get_num_threads`, `openblas_get_num_procs`,
    `openblas_get_config`, `openblas_get_corename` and
    `openblas_get_parallel`.
- `tests/libopenblas_tests.nv` — nine tests over the signatures. Every
  test computes in memory, so the suite reads and writes nothing.

### Double precision, and only double precision

Every routine here is the `d` one. CBLAS names its four precisions with
one letter — `s` for single, `d` for double, `c` and `z` for the two
complex ones — and the single-precision routines take their scalars as
a C `float`.

The novo-lang foreign function interface passes a `Float` as a C
`double`. A C function whose parameter is a `float` reads the low 32
bits of that double, which are the trailing mantissa bits and not a
`float` of the same value; a C function that returns a `float` returns
32 bits where the caller reads 64. Neither direction can be declared.
So `cblas_sgemm`, `cblas_sdot` and the rest of the single-precision
half are absent, and the four groups above are the double-precision
surface in full.

An array of single-precision numbers is not affected: it crosses as an
address like any other buffer. It is the scalars — `alpha`, `beta`, the
rotation's cosine and sine, and every `float` result — that cannot
cross.

### CBLAS, not LAPACK and not the Fortran interface

`libopenblas.so` exports three interfaces: CBLAS, which is C and takes
its arguments by value; the Fortran BLAS, whose symbols end in an
underscore and take every argument including the lengths by pointer;
and LAPACK, with the LAPACKE C interface over it. This package is
CBLAS. The other two are named as missing below.

### The wrapped library

`wraps = "libopenblas"`. Another BLAS — the reference netlib build,
BLIS, Intel's MKL — exports the same CBLAS symbols, so the declarations
resolve against it as well; the six `openblas_` entry points are the
ones that do not.

### Not a `0.0.x` interface release

An interface release is the shape whose every `pub fn` body is a
`todo()`. Every `pub fn` here is an `@ffi` declaration with no body, so
`novo pkg publish` reads the package as a release with bodies and
refuses a `0.0.x` version for it. The first release of a bindings
package is therefore `0.1.0`.

### Unverified

OpenBLAS is not installed on the staging machine, so the suite has
never linked: `novo test` stops at `/usr/bin/ld: cannot find
-lopenblas`. The declarations were checked against the CBLAS interface
as OpenBLAS declares it in `cblas.h`, and `novo pkg build` type-checks
them, which is the whole of what has been measured. Treat the package
as unmeasured until someone runs the suite against a real library.

### Named as missing

**The single-precision half.** `cblas_sgemm`, `cblas_sdot`,
`cblas_snrm2`, `cblas_sasum`, `cblas_saxpy` and every other `s`
routine. The reason is the `float` scalar above, and it applies to all
of them at once.

**The complex halves.** The `c` and `z` routines take their scalars as
a pointer to an interleaved pair rather than by value, so the interface
is expressible; they are left out of the first release because there is
no consumer and no novo-lang type for a complex number.

**LAPACK and LAPACKE.** `LAPACKE_dgesv`, `LAPACKE_dgeqrf` and the rest
of the factorisations live in the same shared object. They are a
separate interface with its own header and its own argument
conventions, and they belong in a package of their own.

**The Fortran interface.** `dgemm_`, `ddot_` and their neighbours take
every argument by pointer, including the dimensions and the transpose
flags, which are single characters. CBLAS is the C interface to the
same kernels.

**`cblas_xerbla`.** The error handler is variadic and is called by the
library rather than by a program. CBLAS reports a bad argument through
it and then, in OpenBLAS, ends the process.

**The banded and packed storage routines.** `cblas_dgbmv`, `cblas_dspmv`,
`cblas_dtpmv` and the rest take a matrix in a storage layout of its own.
They are left out of the first release; the full and triangular layouts
are here.
