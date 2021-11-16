use rand::distributions::Distribution;
use statrs::distribution::Normal;
use rayon::prelude::*;

mod generated {
    type fty = f64;
    type ity = i64;
    type uty = u64;
    
    use std::f64::consts::LOG2_10;
    use std::f64::consts::LOG2_E;
    use std::f64::consts::PI;
    
    fn select(a: bool, b: fty, c: fty) -> fty {
        if a {
            b
        } else {
            c
        }
    }
    
    fn iabs(i: ity) -> ity {
        i.abs()
    }
    
    fn f(f: fty) -> fty {
        f
    }
    
    fn from_bits(u: uty) -> fty {
        fty::from_bits(u)
    }
    
    fn to_bits(f: fty) -> uty {
        fty::to_bits(f)
    }

    fn exp2(arg: fty) -> fty {
        let r: fty = arg.round();
        let mul: fty = fty::from_bits(
            (r.mul_add(0x0010000000000000_u64 as fty, 0x3ff0000000000000_u64 as fty)) as uty,
        );
        let x: fty = arg - r;
        (from_bits(4549839347750377909u64))
            .mul_add(x, from_bits(4563827094295188139u64))
            .mul_add(x, from_bits(4576698039041613846u64))
            .mul_add(x, from_bits(4588159642448921967u64))
            .mul_add(x, from_bits(4597823092488205992u64))
            .mul_add(x, from_bits(4604418534717280147u64))
            .mul_add(x, from_bits(4607182418800017408u64))
            * mul
    }

    fn exp(arg: fty) -> fty {
        exp2(arg * LOG2_E)
    }

    fn negate_on_odd(x: fty, y: fty) -> fty {
        let sign_bit: uty = (((x as ity) & 1) << 63i32) as uty;
        fty::from_bits(sign_bit ^ y.to_bits())
    }
    fn recip_approx(x: fty) -> fty {
        let y: fty = fty::from_bits(
            ((x.abs().to_bits() as fty).mul_add(-1.0, 0x3ff0000000000000_u64 as fty * 2.0)) as uty,
        );
        (y - 0.08).copysign(x)
    }
    fn sqrt_approx(x: fty) -> fty {
        let y: fty = fty::from_bits(
            ((x.abs().to_bits() as fty).mul_add(0.5, 0x3ff0000000000000_u64 as fty * 0.5)) as uty,
        );
        y - 0.08
    }
    fn cbrt_approx(x: fty) -> fty {
        let y: fty = fty::from_bits(
            ((x.abs().to_bits() as fty).mul_add(1.0 / 3.0, 0x3ff0000000000000_u64 as fty * 2.0 / 3.0))
                as uty,
        );
        (y - 0.08).copysign(x)
    }

    pub fn qnorm(arg: fty) -> fty {
        let x: fty = arg - 0.5;
        let recip: fty = 1.0 / (x * x - 0.25);
        let y: fty = (from_bits(4730221388428958202u64))
            .mul_add(x * x, -from_bits(4731626383781768040u64))
            .mul_add(x * x, from_bits(4727627778628654481u64))
            .mul_add(x * x, -from_bits(4720012863723153492u64))
            .mul_add(x * x, from_bits(4708869911609092829u64))
            .mul_add(x * x, -from_bits(4695087533321972728u64))
            .mul_add(x * x, from_bits(4678670384600451257u64))
            .mul_add(x * x, -from_bits(4658680898319303328u64))
            .mul_add(x * x, from_bits(4635605149421499302u64))
            .mul_add(x * x, from_bits(4578476110820645018u64))
            .mul_add(x * x, from_bits(4611041379213747643u64))
            .mul_add(x * x, -from_bits(4603819697584151827u64))
            * x;
        y * recip
    }

    /// See https://xorshift.di.unimi.it/splitmix64.c
    pub fn runif(index: usize) -> f64 {
        let mut z = (index + 1) as u64 * 0x9e3779b97f4a7c15;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
        z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
        z = z ^ (z >> 31);
        from_bits((z >> 2) | 0x3ff0000000000000_u64) - 1.0
    }

    pub fn rnorm(index: usize) -> f64 {
        qnorm(runif(index) * 0.999 + 0.0005)
        // qnorm(runif(index))
    }
}

struct BlackScholes {
    underlying_price : f64,
	strike_price : f64,
	time_to_maturity : f64,
	risk_free_interest_rate : f64,
	volatility_of_the_underlying : f64,
	number_of_iterations : usize,
}

enum OptionType {
    Call,
    Put,
}

impl BlackScholes {
    fn naive_stock_price_function(&self) -> impl Fn(f64) -> f64 {
        let sigma = self.volatility_of_the_underlying;
        let rf = self.risk_free_interest_rate;
        let t = self.time_to_maturity;
        let underlying_price = self.underlying_price;
        let s2 = 0.5 * f64::powi(sigma, 2);
        let sqrt_t = f64::sqrt(t);

        move |rand| underlying_price * f64::exp(t * (rf - s2) + sigma * sqrt_t * rand)
    }

    fn improved_stock_price_function(&self) -> impl Fn(f64) -> f64 {
        let sigma = self.volatility_of_the_underlying;
        let rf = self.risk_free_interest_rate;
        let t = self.time_to_maturity;
        let underlying_price = self.underlying_price;
        let s2 = 0.5 * f64::powi(sigma, 2);
        let sqrt_t = f64::sqrt(t);

        // Reduce the function to its bare bones.
        let k1 = f64::ln(underlying_price) + t * (rf - s2);
        let k2 = sigma * sqrt_t;
        move |rand| f64::exp(k1 + k2 * rand)
    }

    // Naive call option price calculation.
    fn naive_option_price(&self, ot: OptionType) -> f64 {
        let mut r = rand::thread_rng();
        let n = Normal::new(0.0, 1.0).unwrap();

        let stock_price = self.naive_stock_price_function();

        use OptionType::*;
        let average = match ot {
            Call => (0..self.number_of_iterations).map(|_| {
                f64::max(stock_price(n.sample(&mut r)) - self.strike_price, 0.0)
            }).sum::<f64>() / (self.number_of_iterations as f64),
            Put => (0..self.number_of_iterations).map(|_| {
                f64::max(self.strike_price - stock_price(n.sample(&mut r)), 0.0)
            }).sum::<f64>() / (self.number_of_iterations as f64),
        };

        let rf = self.risk_free_interest_rate;
        let t = self.time_to_maturity;
        f64::exp(-1.0 * rf * t) * average
    }

    // improve pricing function
    fn improved_option_price(&self, ot: OptionType) -> f64 {
        let mut r = rand::thread_rng();
        let n = Normal::new(0.0, 1.0).unwrap();

        let stock_price = self.improved_stock_price_function();

        use OptionType::*;
        let average = match ot {
            Call => (0..self.number_of_iterations).map(|_| {
                f64::max(stock_price(n.sample(&mut r)) - self.strike_price, 0.0)
            }).sum::<f64>() / (self.number_of_iterations as f64),
            Put => (0..self.number_of_iterations).map(|_| {
                f64::max(self.strike_price - stock_price(n.sample(&mut r)), 0.0)
            }).sum::<f64>() / (self.number_of_iterations as f64),
        };

        let rf = self.risk_free_interest_rate;
        let t = self.time_to_maturity;
        f64::exp(-1.0 * rf * t) * average
    }

    // Use generated rnorm
    fn rnorm_option_price(&self, ot: OptionType) -> f64 {
        let stock_price = self.improved_stock_price_function();

        use OptionType::*;
        use generated::rnorm;
        let average = match ot {
            Call => (0..self.number_of_iterations).map(|i| {
                f64::max(stock_price(rnorm(i)) - self.strike_price, 0.0)
            }).sum::<f64>() / (self.number_of_iterations as f64),
            Put => (0..self.number_of_iterations).map(|i| {
                f64::max(self.strike_price - stock_price(rnorm(i)), 0.0)
            }).sum::<f64>() / (self.number_of_iterations as f64),
        };

        let rf = self.risk_free_interest_rate;
        let t = self.time_to_maturity;
        f64::exp(-1.0 * rf * t) * average
    }
}

fn main() {
    // println!("{}", (0..1000000).map(|i| generated::runif(i)).sum::<f64>());
    // return;


    let bs = BlackScholes {
        underlying_price : 20.0,
        strike_price : 21.0,
        time_to_maturity : 4.0 / 12.0,
        risk_free_interest_rate : 0.1,
        volatility_of_the_underlying : 0.3,
        number_of_iterations : 100000000,
    };

    // {
    //     let mut f = std::fs::File::create("/tmp/1.csv").unwrap();
    //     use std::io::Write;
    //     writeln!(f, "x,y").unwrap();
    //     for i in 0..100000 {
    //         let x = i as f64 * (1.0/100.0);
    //         // writeln!(f, "{:12.8},{:24.20}", x, generated::qnorm(x)).unwrap();
    //         writeln!(f, "{}, {}", i, generated::rnorm(i)).unwrap();
    //     }
    
    // }

    let step = 2;

    match step {
        0 => {
            println!("call option price = {}", bs.naive_option_price(OptionType::Call));
            println!("put option price  = {}", bs.naive_option_price(OptionType::Put));
        }
        1 => {
            println!("call option price = {}", bs.improved_option_price(OptionType::Call));
            println!("put option price  = {}", bs.improved_option_price(OptionType::Put));
        }
        2 => {
            println!("call option price = {}", bs.rnorm_option_price(OptionType::Call));
            println!("put option price  = {}", bs.rnorm_option_price(OptionType::Put));
        }
        _ => unreachable!(),
    }
}
