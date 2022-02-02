## Results (updated 1 Feb 2022)

These results were collected on a standalone Linux workstation with 32 GB
of memory and a 12-core AMD Threadripper 2920X CPU.

### NAG

* NAG Fortran 7.1.7103, MPICH 3.3.2
* Coarray:
  ```
  nagfor -O3 -coarray -f2018 coarray_collectives.f90 index_map_type.f90 main.f90
  export NAGFORTRAN_NUM_IMAGES=12
  ./a.out <test arguments>
  ```
* MPI:
  ```
  mpifort -O3 index_map_type.f90 main.f90
  mpirun -np 12 -bind-to rr ./a.out <test arguments>
  ```
* Times in msec

NAG   |  MPI |  v1 |  v2 | v3 |  v4 |
------|------|-----|-----|----|-----|
B0-12 |  0.0075 | 0.020 | 0.417 | 0.024 | 0.345 |
B1-12 |  0.010  | 0.022 | 0.420 | 0.028 | 0.333 |
B2-12 |  0.019  | 0.032 | 0.424 | 0.037 | 0.332 |
B3-12 |  0.044  | 0.051 | 0.459 | 0.063 | 0.336 |
B4-12 |  0.085  | 0.129 | 0.524 | 0.098 | 0.385 |


### GNU Fortran / OpenCoarrays

* GCC 11.2.0, MPICH 3.3.2, OpenCoarrays 2.9.2-13-g235167d
* Coarray:
  ```
  caf -O3 coarray_collectives.f90 index_map_type.f90 main.f90
  cafrun -n 12 ./a.out <test arguments>
  ```
* MPI:
  ```
  mpifort -O3 index_map_type.f90 main.f90
  mpirun -np 12 -bind-to rr ./a.out <test arguments>
  ```
* Times in msec

GNU   |  MPI   |  v1 |  v2 |  v3 |  v4 |
------|--------|-----|-----|-----|-----|
B0-12 | 0.0078 |  28 | 9.2 |  31 |  10 |
B1-12 | 0.011  |  58 |  17 |  57 |  18 |
B2-12 | 0.020  | 115 |  33 | 115 |  36 |
B3-12 | 0.045  | 281 |  76 | 272 |  88 |
B4-12 | 0.085  | 489 | 140 | 487 | 160 |


### Intel
* Intel Fortran and MPI 2021.5.0
* Coarray:
  ```
  ifort -O3 -coarray coarray_collectives.f90 index_map_type.f90 main.f90
  export FOR_COARRAY_NUM_IMAGES=12
  export I_MPI_PIN_PROCESSOR_LIST=0-11
  export I_MPI_FABRICS=shm
  export I_MPI_DEVICE=shm
  ./a.out <test arguments>
  ```
* MPI:
  ```
  mpifc -O3 index_map_type.f90 main.f90
  export I_MPI_PIN_PROCESSOR_LIST=0-11
  unset I_MPI_FABRICS I_MPI_DEVICE
  mpirun -np 12 ./a.out <test arguments>
  ```
* Times in msec

Intel |  MPI  |   v1  |  v2   |  v3 |  v4 |
------|-------|-------|-------|-----|-----|
B0-12 | 0.071 |    24 | 0.127 |  17 |   8 |
B1-12 | 0.071 |    79 | 0.151 |  32 |  16 |
B2-12 | 0.081 |   419 | 0.179 |  60 |  29 |
B3-12 | 0.102 |  3100 | 0.244 | 139 |  69 |
B4-12 | 0.169 | 15000 | 0.372 | 257 | 125 |
