# black-scholes-merton

Use rustup https://rustup.rs/ to install cargo and rustc.

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
