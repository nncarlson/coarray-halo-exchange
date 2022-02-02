## Persistent Coarray Buffer for Version 1

NAG comes very close to MPI performance for version 1 of the coarray
implementation. Are there any modifications that would improve its
performance? Note that each invocation of the `gather_aux` procedure
allocates a local coarray buffer `src[*]` (and deallocates it on return).
The amount of storage associated with the coarray is negligible -- merely
a scalar derived type with a single (unassociated) array pointer.
However there may be significant overhead associated with allocation
of a coarray. So two modifications of version 1 were considered that
replaced the local coarray with a persistent coarray that is allocated
once and reused by each invocation of `gather_aux`:
* Version 1A: The local coarray is replaced by a new allocatable component
  `%src[:]` of the `index_map` type. Source in directory `coarray1a`.
* Version 1B: The local coarray is replaced by a module variable.
  Source in directory `coarray1b`.


 time (ms)  |    MPI  |   v1  |  v1A   |  v1B   |
------|---------|-------|--------|--------|
B0-12 |  0.0075 | 0.020 | 0.0050 | 0.0053 |
B1-12 |  0.010  | 0.022 | 0.0077 | 0.0078 |
B2-12 |  0.019  | 0.032 | 0.016  | 0.016  |
B3-12 |  0.044  | 0.051 | 0.035  | 0.037  |
B4-12 |  0.085  | 0.129 | 0.120  | 0.126  |

* There is no real difference between 1A and 1B; not surprising.
* The modified versions are now faster than the original and MPI
  for all but the largest problem.
* The B4 problem is very large with a lot of communication. All
  coarray versions performed about the same. Perhaps some cache effects
  being seen here that are associated with NAG's shared memory coarray
  communication implementation.

### Design Implications

There are significant design downsides to the 1A and 1B approaches.

* For 1B it just feels wrong to have multiple instances of `index_map`
  sharing the same coarray buffer. Perhaps it is okay, but it may be
  prone to "threadsafe" questions.
* For 1B first use can allocate, but who deallocates?
* For 1A, the passed object for the `gather` and `gather_aux` procedures
  must be changed from `intent(in)` to `intent(inout)` because the coarray
  component is being modified. Not terrible, and is something that is always
  accommodated when an otherwise `intent(in)` passed object carries along
  temporary workspace.
* The biggest problem with 1A is that the passed object argument to `init`
  can no longer be `intent(out)`, which is fundamental to the semantics of
  an `init` procedure that is regarded as instantiating an object. This is
  due to a restriction on dummy arguments that have a coarray component --
  apparently the standard writers didn't want allocatable coarray components
  to be implicitly deallocated on entry.
