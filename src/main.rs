use rand::distributions::Distribution;
use statrs::distribution::Normal;

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
    // Naive call option price calculation.
    fn option_price(&self, ot: OptionType) -> f64 {
        let sigma = self.volatility_of_the_underlying;
        let rf = self.risk_free_interest_rate;
        let t = self.time_to_maturity;
        let up = self.underlying_price;
        let s2 = 0.5 * f64::powi(sigma, 2);
        let root_t = f64::sqrt(t);

        let mut r = rand::thread_rng();
        let n = Normal::new(0.0, 1.0).unwrap();

        let stock_price = |rand| up * f64::exp(t * (rf - s2) + sigma * root_t * rand);

        use OptionType::*;
        let average = match ot {
            Call => (0..self.number_of_iterations).map(|_| {
                f64::max(stock_price(n.sample(&mut r)) - self.strike_price, 0.0)
            }).sum::<f64>() / (self.number_of_iterations as f64),
            Put => (0..self.number_of_iterations).map(|_| {
                f64::max(self.strike_price - stock_price(n.sample(&mut r)), 0.0)
            }).sum::<f64>() / (self.number_of_iterations as f64),
        };

        f64::exp(-1.0 * rf * t) * average
    }
}

fn main() {
    let bs = BlackScholes {
        underlying_price : 20.0,
        strike_price : 21.0,
        time_to_maturity : 4.0 / 12.0,
        risk_free_interest_rate : 0.1,
        volatility_of_the_underlying : 0.3,
        number_of_iterations : 100000000,
    };

    println!("call option price = {}", bs.option_price(OptionType::Call));
    println!("put option price  = {}", bs.option_price(OptionType::Put));
}
