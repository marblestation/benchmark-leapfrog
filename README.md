# Simple N-Body with LeapFrog integrator

Implementation in C, Fortran, Go and Rust of a very simple N-Body simulator with 3 particles using a LeapFrog integrator. Presented in [What can the programming language Rust do for astrophysics?](https://arxiv.org/abs/1702.02951), to appear in the Proceedings of the IAU Symposium 325 on Astroinformatics.


## Compilation & execution

From each directory:

```
make clean
make
time target/optimized/leapfrog 
```
