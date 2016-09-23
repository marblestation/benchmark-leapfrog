
all: target/not_optimized/leapfrog target/optimized/leapfrog

target/not_optimized/leapfrog: leapfrog.go target/not_optimized/
	go build -gcflags '-N -l' -o target/not_optimized/leapfrog leapfrog.go
	echo "run: time target/not_optimized/leapfrog"

target/optimized/leapfrog: leapfrog.go target/optimized/
	go build -o target/optimized/leapfrog leapfrog.go
	echo "run: time target/optimized/leapfrog"

target/optimized/:
	mkdir -p target/optimized/

target/not_optimized/:
	mkdir -p target/not_optimized/

clean:
	rm -rf target/not_optimized/leapfrog target/optimized/leapfrog

