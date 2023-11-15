use std::fs::File;
use std::io::Write;
use std::path::Path;

fn main() {
    let out_dir = std::env::var("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("primes.rs");
    let mut file = File::create(dest_path).unwrap();

    let mut primes = vec![];
    let mut primes_until_2k = vec![0];

    'outer: for i in 2u32..=1 << 16 {
        if i.is_power_of_two() {
            primes_until_2k.push(primes.len());
        }
        for p in &primes {
            if i % p == 0 {
                continue 'outer;
            }
        }
        primes.push(i);
    }

    writeln!(
        file,
        "/// First 2<sup>16</sup> prime numbers in a const-sized `[u64; {}]`, generated in \
		 `build.rs`. See more information in crate-level documentation.",
        primes.len()
    )
    .unwrap();
    writeln!(file, "pub static PRIMES: [u64; {}] = [", primes.len()).unwrap();
    for prime in primes {
        writeln!(file, "    {},", prime).unwrap();
    }
    writeln!(file, "];").unwrap();

    writeln!(
        file,
        "/// Mapping from `k` to `π(2<sup>k</sup>)` const-sized `[usize; {}]`, generated in \
		 `build.rs`. See more information in crate-level documentation.",
        primes_until_2k.len()
    )
    .unwrap();
    writeln!(
        file,
        "pub static PRIMES_UNTIL_2K: [usize; {}] = [",
        primes_until_2k.len()
    )
    .unwrap();
    for p in primes_until_2k {
        writeln!(file, "    {},", p).unwrap();
    }
    writeln!(file, "];").unwrap();
}
