# black-scholes-merton

## Initial investigation

This shows how a naive speed improvement can be achieved using a compiled
language:

Use rustup https://rustup.rs/ to install cargo and rustc.

```
$ time python python/bs.py
Call option price with Monte Carlo approach:  1.2405654541052522
Put option price with Monte Carlo approach:  1.5519545170110747

real    0m13.801s
user    0m12.329s
sys     0m2.235s

$ time cargo run --release
    Finished release [optimized] target(s) in 0.01s
     Running `target/release/black-scholes-montecarlo`
call option price = 1.240817518108004
put option price  = 1.5523784672155365

real    0m3.230s
user    0m3.217s
sys     0m0.013s
```

## How to improve this.

### Step 1 is to simplify the function.

The function we evaluate can be reduced to an exponential
of a linear function:

```
        let k1 = f64::ln(underlying_price) + t * (rf - s2);
        let k2 = sigma * sqrt_t;
        move |rand| f64::exp(k1 + k2 * rand)
```

```
$ time target/release/black-scholes-montecarlo 
call option price = 1.2409853902861983
put option price  = 1.5525796331857264

real    0m2.992s
user    0m2.980s
sys     0m0.012s
```

The time is roughly the same (within the noise) showing that
the random number generator is the main source of pain.


### Step 2 is to replace the random number generator with a hash.

Replacing the RNG with a hash enables us to generate 