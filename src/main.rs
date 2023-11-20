use std::convert::identity;
use std::fs::File;
use std::hint::black_box;
use std::io::{stdout, Write};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use e16_primes::*;
use numext_fixed_uint::{u2048 as int, U2048};
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
/// Random int between `a`, `b` with `a < b`
fn random(a: &Int, b: &Int) -> Int {
    let mut out = Int::zero();
    let mut rng = thread_rng();
    let ba = (b - a).0;
    // better branch prediction, handwritten asm?
    let mut i = 0;
    for (j, ba) in ba.iter().enumerate() {
        if ba != &0 {
            i = j;
        }
    }
    out.0[i] = rng.gen_range(0..ba[i]);
    rng.fill(&mut out.0[..i]);
    add(a, &out)
}
/// Checks whether a prime less than `b` divides `a`
///
/// - `b` must be less than or equals 2<sup>16</sub> (enough for `k <= 1024` and correct `C_OPT_D`)
fn trial_division(a: &Int, b: u64) -> bool {
    debug_assert!(b <= 1 << 16, "`b` too big, {b}");
    debug_assert!(&Int::from(b) < a, "`b` greater a, {b} >= {a}");
    PRIMES
        .iter() // avoid `.into_iter()` to prevent stack overflow
        .take(PRIMES_UNTIL_2K[usize::BITS as usize - b.leading_zeros() as usize])
        .all(|p| (a % p).0[0] != 0)
}
/// Generates a number between `0.5` and `max`, see Appendix 1
fn generate_relative_size(max: f64) -> f64 {
    debug_assert!((0.5..=1.0).contains(&max));
    // r = 1 + log2(x)
    // r - 1 = log2(x)
    // 2^(r - 1) = x
    let r = thread_rng().gen_range(0.0..1.0 + max.log2());
    2.0f64.powf(r - 1.0)
}
/// Checks whether lemma 1 is true with `F = q`
fn check_lemma1(k: usize, n: &Int, a: &Int, r: Int, q: &Int) -> bool {
    // (1) a^(n - 1) = 1 mod n
    // (2) gcd(a^((n - 1) / q) - 1, n) = 1
    //
    // (a^((n - 1) / q)) ^ q = a^((n - 1) / q * q) = a^(n - 1)
    //
    // n = 2RF + 1 with F = q
    // n = 2Rq + 1
    // (n - 1) / q = 2R
    // (2) gcd(a^(2R) - 1, n) = 1

    // this is e >> 1
    let e = &r;
    let k_ = (k + 1).div_ceil(u64::BITS as usize);
    let mut r_ = Int::zero();
    r_.0[k_] = 1;
    let r_ = r_;
    let y_ = &r_ % n;

    let mut n_ = mod_inverse(n, &r_);
    n_.0[k_..].fill(0);
    let n_ = sub(&r_, &n_);

    let mut x = Int::zero();
    let l = x.0.len();
    x.0[k_..].copy_from_slice(&a.0[..l - k_]);
    x %= n;
    let mut y = y_.clone();
    for i in 0..Int::count_bits() as usize - e.leading_zeros() as usize - 1 + 1 {
        // only least significant bit relevant for odd check, sub 1 since first is 0
        if e.bit(i.wrapping_sub(1)) == Some(true) {
            y = redc(k_, n, &n_, mul(&x, &y));
        }
        x = redc(k_, n, &n_, mul(&x, &x));
    }

    let anq = redc(k_, n, &n_, mul(&x, &y));

    if &(&anq - Int::one()).gcd(&y_) == &y_ {
        return false;
    }

    let e = q;
    let mut x = anq;
    let mut y = y_.clone();
    for i in 0..Int::count_bits() as usize - e.leading_zeros() as usize - 1 {
        // only least significant bit relevant for odd check
        if e.bit(i) == Some(true) {
            y = redc(k_, n, &n_, mul(&x, &y));
        }
        x = redc(k_, n, &n_, mul(&x, &x));
    }

    let anq2 = redc(k_, n, &n_, mul(&x, &y));

    anq2 == y_
}

// constants
const C_OPT_N: u64 = 1;
const C_OPT_D: u64 = 1024 * 1024 / (1 << 16);
const MARGIN: u64 = 16;
const P0: u64 = 2 * MARGIN; // 2 * margin; if less than (2^16)^2 then easy to find prime with table

fn fast_prime(k: u64) -> Option<Int> {
    let n;

    if k <= P0 {
        let _2k = 1 << k;
        let _2k1 = _2k >> 1;
        let mut n64 = 4; // not prime
        let mut rng = thread_rng();
        while PRIMES
            .iter()
            .take(PRIMES_UNTIL_2K[(k as usize - 1) / 2 + 1])
            .any(|p| n64 % p == 0)
        {
            n64 = rng.gen_range(_2k1.._2k);
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
        #[allow(non_snake_case)]
        let I2 = &I << 1;
        n = std::iter::repeat(false)
            // .take(2000) // use this to limit tries, this results in fails
            .par_bridge()
            .into_par_iter()
            .map(|_| {
                let r = random(&I, &I2);
                let mut n = mul(&r, &q2);
                let a = random(&int!("2"), &n);
                n.0[0] += 1; // this can't overflow since rq2 is even (u64::MAX is odd)

                (trial_division(&n, g) && check_lemma1(k as usize, &n, &a, r, &q)).then_some(n)
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

// R = 2^(k * 64)
fn redc(k: usize, n: &Int, n_: &Int, t: Int) -> Int {
    let mut t2 = t.clone();
    // t %= R
    t2.0[k..].fill(0);
    // m = (t * n') % R
    let mut m = mul(&t2, n_);
    m.0[k..].fill(0);

    // t' = (t + m * n) / R
    let mn = mul(&m, n);
    let mut t_ = Int::zero();
    let mut borrow = false;
    for i in 0..k {
        let (a, b) = t.0[i].overflowing_add(mn.0[i]);
        borrow = b || a.overflowing_add(1).1;
    }
    for i in k..t_.0.len() {
        let b;
        (t_.0[i - k], b) = t.0[i].overflowing_add(mn.0[i]);
        let mut b2 = false;
        if borrow {
            (t_.0[i - k], b2) = t_.0[i - k].overflowing_add(1);
        }
        borrow = b || b2;
    }
    while &t_ > n {
        t_ = sub(&t_, n);
    }
    t_
}
fn extended_gcd(a: &Int, b: &Int) -> (Int, Int, bool) {
    let (mut old_r, mut r) = (a.clone(), b.clone());
    let (mut old_s, mut s, mut sign_s, mut sign_old_s) = (Int::one(), Int::zero(), false, false);

    while !r.is_zero() {
        let quotient = &old_r / &r;
        let qr = mul(&quotient, &r);
        let qs = mul(&quotient, &s);

        (old_r, r) = (r, &old_r - qr);
        let s_ = (s.clone(), sign_s);
        (s, sign_s) = match (sign_s, sign_old_s) {
            (false, false) if old_s >= qs => (sub(&old_s, &qs), false),
            (false, false) => (sub(&qs, &old_s), true),
            (false, true) => (add(&old_s, &qs), true),
            (true, false) => (add(&old_s, &qs), false),
            (true, true) if qs >= old_s => (sub(&qs, &old_s), false),
            (true, true) => (sub(&qs, &old_s), true),
        };
        (old_s, sign_old_s) = s_;
    }

    (old_r, old_s, sign_old_s)
}

fn mod_inverse(a: &Int, n: &Int) -> Int {
    let (gcd, x, s) = extended_gcd(a, n);
    debug_assert_eq!(gcd, Int::one());
    if s {
        sub(&n, &x)
    } else {
        x
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mod_inverse() {
        let a: u32 = 67;
        let n: u32 = 101;
        let e: u32 = 98;
        assert_eq!(Int::from(e), mod_inverse(&Int::from(a), &Int::from(n)));
    }
}
