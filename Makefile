# CC=g++
# CFLAGS=-mcmodel=medium -fopenmp -O -march=native -fopt-info-vec-optimized -mavx -mavx2
# # compiler may crash when static array too large,
# # add `-mcmodel=medium` in this case.

# all:
# 	$(CC) -o gemm hw_baseline.cpp $(CFLAGS)

# .PHONY: run
# run: all
# 	./gemm

# .PHONY: clean
# clean:
# 	rm -rf *.o gemm



CC=mpiicc
CFLAGS=-mcmodel=medium -fopenmp -Ofast -march=native -mavx -mavx2
# compiler may crash when static array too large,
# add `-mcmodel=medium` in this case.

all:
	$(CC) -o gemm_mpi mpi_async.cpp $(CFLAGS)

.PHONY: run
run: all
	mpirun -ppn 1 ./gemm_mpi



# CC=mpicc
# CFLAGS=-mcmodel=medium -fopenmp -Ofast -march=native -mavx -mavx2
# # compiler may crash when static array too large,
# # add `-mcmodel=medium` in this case.

# all:
# 	$(CC) -o gemm_mpi mpi_async_gcc.cpp $(CFLAGS)

# .PHONY: run
# run: all
# 	mpirun -ppn 1 ./gemm_mpi