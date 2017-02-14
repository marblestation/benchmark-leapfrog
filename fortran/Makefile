
all: target/not_optimized/leapfrog target/optimized/leapfrog

target/not_optimized/leapfrog: leapfrog.f90 target/not_optimized/
	gfortran -O0 -finit-real=nan -fdefault-real-8 leapfrog.f90 -o target/not_optimized/leapfrog
	echo "run: time target/not_optimized/leapfrog"

target/optimized/leapfrog: leapfrog.f90 target/optimized/
	gfortran -O3 -finit-real=nan -fdefault-real-8 leapfrog.f90 -o target/optimized/leapfrog
	echo "run: time target/optimized/leapfrog"

target/optimized/:
	mkdir -p target/optimized/

target/not_optimized/:
	mkdir -p target/not_optimized/

clean:
	rm -rf target/not_optimized/leapfrog target/optimized/leapfrog

