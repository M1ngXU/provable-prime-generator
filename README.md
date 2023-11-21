# Provable primes generator

A Rust implementation to generate provable primes using [Maurer's Method](https://crypto.ethz.ch/publications/files/Maurer95a.pdf).

## Implementation

This implementation is optimized for 1024-bit primes, using the "Simplified Version" (3.4). The branch `montgomery-redc` uses [Montgomery Reduction](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication) to exponentiate modulo `N` faster; but with only 1024 bit this method is slightly slower.

<i>Note: this branch is slightly outdated</i>

## Execute
After cloning the repository, running `cargo run --release` measures the average execution time to generate a 1024-bit prime and outputs the duration of each generation in `stats.tsv`. If the maximum running time of each recursion is set, then a generation might fail (i.e. no prime was found after a certain amount of iterations) and the first column of `stats.tsv` contains an `f` (for fail). On success the first column is an `s`.

Adding the feature `ignore-overflow` silently wraps overflows; since this should not happen, this gives some speedup compared to the default overflow-checking.

By default this uses multithreading with `rayon` to speed up the search for prime numbers.
