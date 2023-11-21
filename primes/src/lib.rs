//! First 2<sup>16</sup> prime numbers in a const-sized [ `u64` ] array,
//! generated in build.rs. Also includes a mapping from `k` to `π(2<sup>k</sup>)`.

include!(concat!(env!("OUT_DIR"), "/primes.rs"));
