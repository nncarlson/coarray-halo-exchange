### Intel Issues

Here are some specific issues with the Intel compiler that I encountered

* By default the coarray tests were needlessly hammering the network
  (loopback), which resulted in *even larger* times than now reported.
  This was solved by setting the environment variables `I_MPI_FABRICS`
  and `I_MPI_DEVICE` to `shm`; see this
  [Intel page](https://www.intel.com/content/www/us/en/develop/documentation/fortran-compiler-oneapi-dev-guide-and-reference/top/optimization-and-programming-guide/coarrays-1/using-coarrays.html).

* The mpi tests will fail to run if either `I_MPI_DEVICE` or `I_MPI_FABRICS`
  are defined. An error message to that effect is written.

* The mpi tests hammer on the network (loopback) just like the coarray tests
  but here I have no solution.
