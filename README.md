# libopenblas-sys

BLAS, the Basic Linear Algebra Subprograms, is a specification of the
operations a program performs on vectors and matrices: a dot product, a
matrix times a vector, a matrix times a matrix. It is specified by
[the BLAS Technical Forum standard](https://www.netlib.org/blas/blast-forum/),
and its C form is
[CBLAS](https://www.netlib.org/blas/#_cblas). OpenBLAS is an
implementation of it whose inner loops are written in assembly for each
processor family, and it is documented at
[openblas.net](https://www.openblas.net/). This package declares
twenty-seven of that library's entry points to novo-lang, one
declaration each.

**Status: a binding, not a port.** Every function in this package is a
declaration of a function in libopenblas. The package contains no logic
of its own, and it does nothing without the C library installed. The
twenty-seven entry points are the double-precision BLAS routines and
the six queries OpenBLAS adds to them; the section "What is not
included" says what a program still cannot do with them alone.

**Unverified.** OpenBLAS is not installed on the machine where this
package was written, so the test suite has never been linked. See the
"Tests" section.

## What it is

A **vector** is a run of numbers in memory, described by three things:
where it starts, how many elements it has, and the stride from one
element to the next. A stride of 1 walks a contiguous run. A stride of
2 walks every second element, which is how a program reads one column
out of a matrix stored by rows.

A **matrix** is a run of numbers with a shape. It is stored either by
rows, so that the elements of one row are next to each other, or by
columns. The **leading dimension** is how far apart the starts of two
neighbouring rows are, counted in elements. It is at least the row
length and it may be more, which is what lets a program hand a
sub-block of a larger matrix to a routine without copying it.

The routines are grouped into three **levels**, and the level says how
much work each one does for the data it reads. Level 1 is a vector
against a vector. Level 2 is a matrix against a vector. Level 3 is a
matrix against a matrix, and it is where a tuned implementation is
many times faster than an obvious one, because each element read is
used many times.

Almost every routine computes into an existing result rather than
returning a new one. The shape is `C = alpha * (something) + beta * C`:
`alpha` scales the product, `beta` scales what the destination already
holds, and a `beta` of zero overwrites it. The destination is one of
the arguments, and it is written in place.

A routine's name carries its precision in its first letter. `d` is
double precision, which is the eight-byte IEEE-754 number novo-lang
calls a `Float`. Every routine in this package is a `d` routine. The
single-precision `s` routines are absent, and the reason is in "What is
not included".

## Install

```
novo pkg add libopenblas-sys
```

Adding the package does not install the C library. On Debian and Ubuntu
the library and its header come from the system package
`libopenblas-dev`:

```
sudo apt install libopenblas-dev
```

On macOS the Homebrew formula is `openblas`. On other systems the
library builds from the OpenBLAS source with `make`.

## Example

A two by two matrix multiplied by a vector:

```novo ignore
use libopenblas

fn main() [io, ffi]
    // A = [[1, 2], [3, 4]], stored by rows, so the leading dimension
    // is 2.
    let a = ptr.alloc(4 * 8)
    ptr.write_word(a, float.to_bits(1.0))
    ptr.write_word(a + 8, float.to_bits(2.0))
    ptr.write_word(a + 16, float.to_bits(3.0))
    ptr.write_word(a + 24, float.to_bits(4.0))

    // x = [1, 1], and y is the destination.
    let x = ptr.alloc(2 * 8)
    ptr.write_word(x, float.to_bits(1.0))
    ptr.write_word(x + 8, float.to_bits(1.0))
    let y = ptr.alloc(2 * 8)

    // 101 is row-major storage and 111 leaves A as it stands.
    // A beta of 0.0 overwrites y rather than adding to it.
    libopenblas.cblas_dgemv(101, 111, 2, 2, 1.0, a, 2, x, 1, 0.0, y, 1)
    println("${float.from_bits(ptr.read_word(y))}")
    println("${float.from_bits(ptr.read_word(y + 8))}")

    ptr.free(y)
    ptr.free(x)
    ptr.free(a)
```

The example is fenced as an illustration rather than a compiled block
because `novo doc` compiles the blocks in documentation comments and
not the ones in this file. The same calls are in
`tests/libopenblas_tests.nv`.

## What the package contains

| Module | Contents |
| --- | --- |
| `libopenblas` | Every entry point, in four groups: the level 1 vector routines, the level 2 matrix-vector routines, the level 3 matrix-matrix routines, and the six OpenBLAS runtime queries. |

The four groups and their sizes:

| Group | Entry points | What it does |
| --- | --- | --- |
| Level 1 | 10 | Dot product, the two norms, the position of the largest element, scaling, copying, swapping, accumulation and the plane rotation. |
| Level 2 | 5 | A matrix times a vector, in the general, symmetric and triangular forms, and the outer product. |
| Level 3 | 6 | A matrix times a matrix, in the general, symmetric and triangular forms, and the two rank updates. |
| OpenBLAS runtime | 6 | Sets and reads the thread count, counts the processors, and names the build and the processor kernel. |

## How to choose an entry point

`cblas_dgemm` is the general matrix product, and it is the routine to
reach for unless one of the operands has a shape worth naming.

`cblas_dsymm`, `cblas_dsyrk` and `cblas_dsyr2k` are for a symmetric
operand. They read one triangle and do about half the arithmetic of the
general routine.

`cblas_dtrmm` and `cblas_dtrsm` are for a triangular operand.
`cblas_dtrsm` solves rather than multiplies, and it is how the result
of a factorisation is applied to a right-hand side.

The level 2 routines have the same four shapes for a single vector. Use
them when there is one right-hand side rather than a block of them; use
level 3 when there are several, because one level 3 call is faster than
the same work as a loop of level 2 calls.

`openblas_set_num_threads` belongs in a program that parallelises above
the BLAS. Set it to 1 there, or the two layers of threads multiply.

## The rules a user needs

1. **A vector or a matrix is an address, and an element is eight
   bytes.** `ptr.alloc(n * 8)` reserves `n` elements, and the routines
   read and write them in place.
2. **A `Float` is written into memory through its bits.** There is no
   `ptr.write_f64`. `ptr.write_word(addr, float.to_bits(v))` writes one
   element and `float.from_bits(ptr.read_word(addr))` reads it back.
   `float.to_bits` reinterprets rather than converts, so the round trip
   is exact.
3. **The order, transpose, triangle, diagonal and side arguments are
   numbers**, because the C header spells them as enumerations.

   | Argument | Number | Meaning |
   | --- | --- | --- |
   | `CblasRowMajor` | 101 | the elements of a row are next to each other |
   | `CblasColMajor` | 102 | the elements of a column are next to each other |
   | `CblasNoTrans` | 111 | use the operand as it stands |
   | `CblasTrans` | 112 | transpose the operand first |
   | `CblasConjTrans` | 113 | the complex case; for a real matrix it is 112 |
   | `CblasUpper` | 121 | the upper triangle holds the matrix |
   | `CblasLower` | 122 | the lower triangle holds it |
   | `CblasNonUnit` | 131 | the diagonal is stored |
   | `CblasUnit` | 132 | every diagonal element is 1, whatever is stored |
   | `CblasLeft` | 141 | the named operand multiplies from the left |
   | `CblasRight` | 142 | it multiplies from the right |

4. **The leading dimension is in elements, not bytes**, and for
   row-major storage it is the distance between the starts of two
   rows. For a whole `m` by `n` matrix stored by rows it is `n`.
5. **`cblas_idamax` counts from zero.** The Fortran routine it is
   layered over counts from one. This is the one place where the C
   interface changes an answer rather than an argument order.
6. **A destination is read as well as written unless `beta` is zero.**
   A destination that has just been allocated holds whatever was in
   memory, so pass a `beta` of 0.0 or clear it first. `ptr.alloc` does
   not clear.
7. **The triangle that is not named is neither read nor written.**
   `cblas_dsyrk` writes one triangle of its answer and leaves the other
   one as it found it. A program that needs the whole matrix copies the
   triangle across itself.
8. **A solve does not check for a singular matrix.** `cblas_dtrsv` and
   `cblas_dtrsm` divide by the diagonal. A zero there gives an infinity
   or a not-a-number rather than an error.
9. **An entry point that answers a C `int` answers it in 32 bits.**
   Write `as i32` before comparing the answer with a negative number.
   The four `openblas_get_` entry points that answer a count are the
   ones that need it.
10. **A length, a stride and a leading dimension are `blasint`**, which
    is a C `int` in the ordinary build of OpenBLAS and 64 bits wide in
    a build whose `openblas_get_config` string contains `USE64BITINT`.
    A count is positive, so it crosses unchanged either way.
11. **A bad argument ends the process.** CBLAS reports it through
    `cblas_xerbla`, and OpenBLAS's implementation of that prints a line
    and exits. Check the dimensions before the call; there is no return
    code to test.
12. **The thread count is global.** `openblas_set_num_threads` applies
    to every later call from every thread in the process.

## What is not included

- **The single-precision routines.** `cblas_sgemm`, `cblas_sdot` and
  every other `s` routine take their scalars as a C `float`. The
  novo-lang foreign function interface passes a `Float` as a C
  `double`, so a `float` parameter receives the trailing mantissa bits
  of the double and a `float` result is read as a double. An array of
  single-precision numbers is unaffected — it crosses as an address —
  but the scalars cannot cross, so the routines cannot be declared.
- **The complex routines.** The `c` and `z` routines take their scalars
  as a pointer to an interleaved pair of numbers, so the interface is
  expressible. They are left out of the first release.
- **LAPACK and LAPACKE.** `LAPACKE_dgesv`, `LAPACKE_dgeqrf` and the
  factorisations are in the same shared object, under a separate
  interface with its own header. They belong in a package of their own.
- **The Fortran interface.** `dgemm_`, `ddot_` and their neighbours
  take every argument by pointer, including the dimensions and the
  single-character transpose flags. CBLAS is the C interface to the
  same kernels.
- **`cblas_xerbla`.** It is variadic, and the library calls it rather
  than a program.
- **The banded, packed and Hermitian storage layouts.**
  `cblas_dgbmv`, `cblas_dspmv`, `cblas_dtpmv` and their neighbours take
  a matrix in a storage layout of its own. The full, symmetric and
  triangular layouts are here.

## Related packages

`linalg-nv` is linear algebra written in novo-lang, with no C library.
It is the port, and it is the reference for what an answer should be.
This package is the implementation that gets that answer quickly: on a
matrix product of any size, a tuned BLAS is several times faster than
a straightforward loop, and the gap widens with the size.

Choose `linalg-nv` when the program must build for a microcontroller or
for WebAssembly, when a C toolchain is not wanted, or when the matrices
are small enough that the call overhead dominates. Choose this package
when the program spends its time in matrix arithmetic. `linalg-nv` is
planned and not published yet.

`ndarray-nv` is the array container the two of them operate on.

## Tests

`tests/libopenblas_tests.nv` holds nine tests written against the
signatures. They call the C library, so `novo test` needs OpenBLAS
installed and linkable:

```
novo test tests/libopenblas_tests.nv
```

`novo pkg build` type-checks the declarations and needs nothing
installed.

Every test computes in memory, so the suite reads and writes nothing.
The level 1 tests assert the dot product, the two norms and the
zero-based position of the largest element over the vectors 1 to 4 and
5 to 8, and that a stride of 2 reads half of them. The rotation test
asserts that the rotation computed for the point (3, 4) has length 5 and
leaves nothing on the second axis. The level 2 tests multiply a two by
two matrix by a vector, transposed and not, and assert that the
symmetric form reads only the named triangle. The level 3 tests
multiply by the identity and assert the matrix comes back, and that a
rank update writes one triangle and leaves the other at zero. The two
triangular tests multiply and then solve, and assert that the vector
comes back to what it started as.

**Unverified.** OpenBLAS is not installed on the staging machine, so
the suite has never linked: `novo test` stops at `/usr/bin/ld: cannot
find -lopenblas`. Every assertion above is written from the CBLAS
interface and none of them has been observed to pass.

## Implementation status

| Group | State |
| --- | --- |
| Level 1, double precision | Complete. |
| Level 2, double precision | Complete for the full, symmetric and triangular layouts. |
| Level 3, double precision | Complete for the full, symmetric and triangular layouts. |
| OpenBLAS runtime | Complete. |
| Single precision | Absent. A C `float` scalar cannot cross the foreign function interface. |
| Complex precisions | Absent. Left out of the first release. |
| Banded and packed layouts | Absent. Left out of the first release. |
| LAPACK | Absent, and a separate interface. |
| Fortran interface | Absent, and a separate interface. |

## Licence

Apache-2.0. See [LICENSE](LICENSE).

OpenBLAS itself is distributed under the BSD three-clause licence, and
installing it is the reader's own step.
