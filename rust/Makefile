
all: target/not_optimized/leapfrog target/optimized/leapfrog

target/not_optimized/leapfrog: leapfrog.rs target/not_optimized/
	rustc -C opt-level=0 leapfrog.rs -o target/not_optimized/leapfrog
	echo "run: time target/not_optimized/leapfrog"

target/optimized/leapfrog: leapfrog.rs target/optimized/
	rustc -C opt-level=3 leapfrog.rs -o target/optimized/leapfrog
	echo "run: time target/optimized/leapfrog"

target/optimized/:
	mkdir -p target/optimized/

target/not_optimized/:
	mkdir -p target/not_optimized/

clean:
	rm -rf target/not_optimized/leapfrog target/optimized/leapfrog

