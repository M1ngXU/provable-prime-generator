# provable-prime-generator

A Rust implementation to generate provable primes using [Maurer's Method](https://crypto.ethz.ch/publications/files/Maurer95a.pdf).

## Implementation

This implementation is optimized for 1024-bit primes, using the "Simplified Version" (3.4). The branch `montgomery-redc` uses [Montgomery Reduction](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication) to exponentiate modulo `N` faster; but with only 1024 bit this method is slightly slower.
