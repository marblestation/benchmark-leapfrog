
all: target/not_optimized/leapfrog target/optimized/leapfrog

target/not_optimized/leapfrog: leapfrog.c target/not_optimized/
	gcc -std=c99 -O0 -Wall leapfrog.c -lm -o target/not_optimized/leapfrog
	echo "run: time target/not_optimized/leapfrog"

target/optimized/leapfrog: leapfrog.c target/optimized/
	gcc -std=c99 -O3 -Wall leapfrog.c -lm -o target/optimized/leapfrog
	echo "run: time target/optimized/leapfrog"

target/optimized/:
	mkdir -p target/optimized/

target/not_optimized/:
	mkdir -p target/not_optimized/

clean:
	rm -rf target/not_optimized/leapfrog target/optimized/leapfrog

