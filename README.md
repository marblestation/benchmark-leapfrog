# Simple N-Body with LeapFrog integrator

Implementation in C, Fortran, Go, Julia and Rust of a very simple N-Body simulator with 3 particles using a LeapFrog integrator. Presented in [What can the programming language Rust do for astrophysics?](https://arxiv.org/abs/1702.02951), to appear in the Proceedings of the IAU Symposium 325 on Astroinformatics.

The code has evolved since publication, implementing several suggestion made by the [Hacker News](https://news.ycombinator.com/item?id=13632894) and [Rust subreddit](https://www.reddit.com/r/rust/comments/5trref/what_can_rust_do_for_astrophysics/) communities. The current times on a 1,6 GHz Intel Core i5 machine:

- C: 2m53.504s
- Fortran: 3m16.314s
- Rust: 2m33.082s
- Go: 4m10.233s
- Julia: 3m59.223s

Output positions for the two particles after a one million year simulation:

```
[
    [0.00011259641961937496, 0.00011235312238407964, 0.00000982962683070781],
    [-2.986377308237493, 14588403.613911422, 1276320.0322206458]
]
```

The original article's conclusions are still valid for this simple user case, Rust can be as fast and precise as C or Fortran. Additionally, Rust characteristics ensures that scientific results are not affected by memory management issues and this is a great advantage for reliable scientific computation.

## Compilation & execution

From each directory:

```
make clean
make
time target/optimized/leapfrog 
```
