#![feature(concat_idents, portable_simd)]

use std::convert::identity;
use std::fs::File;
use std::hint::black_box;
use std::io::{stdout, Write};
use std::simd::Simd;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use numext_fixed_uint::U2048;
use primes::*;
use rand::{thread_rng, Rng};
use rayon::iter::{IntoParallelIterator, ParallelBridge, ParallelIterator};

type Int = U2048;

fn main() {
    let mut stats = File::create("stats.tsv").unwrap();

    let last_ctrlc = Arc::new(Mutex::new(Instant::now()));
    ctrlc::set_handler(move || {
        stdout().flush().unwrap();

        let now = Instant::now();
        let mut lock = last_ctrlc.lock().unwrap();
        let last = *lock;
        *lock = now;

        if now - last < Duration::from_secs(2) {
            println!();
            std::process::exit(0);
        }
    })
    .expect("Error setting Ctrl-C handler");

    println!("Set Ctrl-C handler");

    let mut total_elapsed = Duration::ZERO;
    let mut fails_elapsed = Duration::ZERO;
    let mut fails = 0;

    for i in 1.. {
        let start = Instant::now();
        let prime = black_box(fast_prime(1023));
        let e = start.elapsed();
        if prime.is_none() {
            fails += 1;
            fails_elapsed += e;
            writeln!(stats, "f\t{e:?}").unwrap();
        } else {
            total_elapsed += e;
            writeln!(stats, "s\t{e:?}").unwrap();
        }
        print!(
            "Average elapsed: {:.2?}/{:.2?}  (runs: {i}, fails: {fails})                \r",
            total_elapsed / (i - fails).max(1),
            fails_elapsed / fails.max(1)
        );
    }
}

/// Calculates 2<sup>a</sup> wrapping around `Int::MAX`
fn power2(a: u64) -> Int {
    let mut out = Int::zero();
    out.set_bit(a as usize, true);
    out
}
/// Random int between `0` and `a`
fn random(a: &Int) -> Int {
    let mut out = Int::zero();
    let mut rng = thread_rng();
    // better branch prediction, handwritten asm?
    let mut i = 0;
    for (j, ba) in a.0.iter().enumerate() {
        if ba != &0 {
            i = j;
            // break; // this line adds 20% to the runtime?? (also add `.rev()` after `.enumerate()`)
        }
    }
    out.0[i] = rng.gen_range(0..a.0[i]);
    rng.fill(&mut out.0[..i]);
    out
}
/// Checks whether a prime less than `b` divides `a`
///
/// - `b` must be less than or equals 2<sup>P0/2</sub> (this actually only checks for primes less than the last power of 2 of b and less than P0/2)
fn trial_division(a: &Int, b: u64) -> bool {
    debug_assert!(b <= 1 << (P0 / 2), "`b` too big, {b}");
    debug_assert!(&Int::from(b) < a, "`b` greater a, {b} >= {a}");
    PRIMES
        .iter()
        .take(PRIMES_UNTIL_2K[u64::BITS as usize - b.leading_zeros() as usize - 1])
        .all(|p| (a % p).0[0] != 0) // p < u64::MAX so fine
}
/// Generates a number between `0.5` and `max`, see Appendix 1
fn generate_relative_size(max: f64) -> f64 {
    // return 0.5; // sacrifice prime distribution for correct (probably really bad)
    if max <= 0.5 {
        0.5
    } else {
        // r = 1 + log2(x)
        // r - 1 = log2(x)
        // 2^(r - 1) = x
        let r = thread_rng().gen_range(0.0..1.0 + max.log2());
        2.0f64.powf(r - 1.0)
    }
}
/// Calculates a<sup>n + additional</sup> mod m
fn pow_mod_n(mut x: Int, n: &Int, m: &Int, additional: usize) -> Int {
    let mut y = Int::one();
    for i in 0..Int::count_bits() as usize - n.leading_zeros() as usize - 1 + additional {
        // only least significant bit relevant for odd check
        //
        // if underflow then the index is too big and returns none (would be not odd in that case)
        if n.bit(i.wrapping_sub(additional)) == Some(true) {
            y = mul(&x, &y) % m;
        }
        x = mul(&x, &x) % m;
    }

    mul(&x, &y) % m
}

/// Checks whether lemma 1 is true with `F = q`
fn check_lemma1(n: &Int, a: Int, r: Int, q: &Int) -> bool {
    // (1) a^(n - 1) = 1 mod n
    // (2) gcd(a^((n - 1) / q) - 1, n) = 1
    //
    // (a^((n - 1) / q)) ^ q = a^((n - 1) / q * q) = a^(n - 1)
    //
    // n = 2RF + 1 with F = q
    // n = 2Rq + 1
    // (n - 1) / q = 2R
    // (2) gcd(a^(2R) - 1, n) = 1

    let anq = pow_mod_n(a, &r, n, 1);
    let anq1 = sub(&anq, &Int::one());

    Simd::from(pow_mod_n(anq, q, n, 0).0) == Simd::from(Int::one().0)
        && Simd::from(anq1.gcd(n).0) == Simd::from(Int::one().0)
}

// constants
const C_OPT_N: u64 = 1;
const C_OPT_D: u64 = 200;
const MARGIN: u64 = 20;

static PRIMES: &[u64] = &PRIMES_10;
const P0: u64 = 2 * 10;

fn fast_prime(k: u64) -> Option<Int> {
    let n;

    if k <= P0 {
        let _2k = 1 << k;
        let _2k1 = _2k >> 1;
        let mut n64; // not prime
        let mut rng = thread_rng();
        loop {
            n64 = rng.gen_range(_2k1.._2k);
            if PRIMES
                .iter()
                .take(PRIMES_UNTIL_2K[(k as usize - 1) / 2 + 1])
                .any(|p| n64 % p == 0)
            {
                break;
            }
        }
        n = Int::from(n64);
    } else {
        // k * relative_size < k - margin
        // margin < k - k * relative_size
        // margin < k * (1 - relative_size)
        // margin < k * (relative_size_d - relative_size_n) / relative_size_d
        // margin * relative_size_d < k * (relative_size_d - relative_size_n)
        //
        // relative_size < (k - margin) / k
        // relative_size < 1 - margin / k
        //
        // edge case:
        //    relative_size = 0.5 (minimum):
        //        => margin < k - 0.5k = 0.5k => margin must be less than 0.5k
        //    if margin is (almost) 0.5k, then it's hard to find a relative_size
        let ratio = MARGIN as f64 / k as f64;
        let relative_size = generate_relative_size(1.0 - ratio);

        let q = fast_prime((relative_size * k as f64).trunc() as u64)?;
        let q2 = &q << 1;
        #[allow(non_snake_case)]
        let I = power2(k - 1) / &q;
        let g = C_OPT_N * k * k / C_OPT_D;
        n = std::iter::repeat(false)
            // .take(2000) // use this to limit tries, this results in fails
            .par_bridge()
            .into_par_iter()
            .map(|_| {
                let r = add(&I, &random(&I));
                let mut n = mul(&r, &q2);
                let mut a = random(&n);
                a.0[0] |= 0b10; // make sure that `a` is at least 2, so a \in [2, n)
                n.0[0] += 1; // this can't overflow since rq2 is even (u64::MAX is odd)

                (trial_division(&n, g) && check_lemma1(&n, a, r, &q)).then_some(n)
            })
            .find_map_any(identity)?;
        // .find_map(identity)?; // use this for single threaded
    }
    Some(n)
}

#[inline(always)]
fn add(a: &Int, b: &Int) -> Int {
    #[cfg(not(feature = "ignore-overflow"))]
    return a + b;
    #[cfg(feature = "ignore-overflow")]
    return a.overflowing_add(b).0;
}
#[inline(always)]
fn sub(a: &Int, b: &Int) -> Int {
    #[cfg(not(feature = "ignore-overflow"))]
    return a - b;
    #[cfg(feature = "ignore-overflow")]
    return a.overflowing_sub(b).0;
}
#[inline(always)]
fn mul(a: &Int, b: &Int) -> Int {
    #[cfg(not(feature = "ignore-overflow"))]
    return a * b;
    #[cfg(feature = "ignore-overflow")]
    return a.overflowing_mul(b).0;
}
