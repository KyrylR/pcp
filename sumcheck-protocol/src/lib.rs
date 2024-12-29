extern crate core;

use crate::fp::Fp;
use crate::mpolynomial::MPolynomial;
use crate::upolynomial::UPolynomial;

pub mod fp;
pub mod mpolynomial;
pub mod upolynomial;

pub struct SumcheckProtocol {
    polynomial: MPolynomial,
    interaction_completed: bool,
    randomness: Vec<Fp>,
    statement: Fp,
    step: usize,
    previous_step: Option<UPolynomial>
}

impl SumcheckProtocol {
    pub fn new(polynomial: MPolynomial) -> Self {
        Self {
            polynomial,
            interaction_completed: false,
            randomness: vec![],
            statement: Fp(0),
            step: 0,
            previous_step: None
        }
    }

    pub fn prove(&mut self) -> Option<UPolynomial> {
        if self.step > self.polynomial.number_of_vars() || self.interaction_completed {
            return None;
        }

        if self.statement == Fp(0) {
            self.statement = self.polynomial.sum_over_hyper_cube(None);
        }

        Some(self.polynomial.fix_var_over_hyper_cube(Some(&self.randomness)))
    }

    pub fn verify(&mut self, rec_polynomial: Option<UPolynomial>) {
        assert!(rec_polynomial.is_some());

        if self.interaction_completed {
            return;
        }

        assert!(rec_polynomial.clone().unwrap().degree() <= self.polynomial.degree_ind());

        // FIXME: consider bigger degrees
        let res1 = rec_polynomial.clone().unwrap().eval(&*vec![Fp(0)]);
        let res2 = rec_polynomial.clone().unwrap().eval(&*vec![Fp(1)]);

        if self.previous_step.is_none() {
            self.previous_step = rec_polynomial.clone();
            assert_eq!(res1 + res2, self.statement);
        } else {
            assert!(self.randomness.len() > 0);
            let latest_randomness = self.randomness.get(self.randomness.len() - 1).unwrap();
            assert_eq!(res1 + res2, self.previous_step.clone().unwrap().eval(&*vec![*latest_randomness]));
            self.previous_step = rec_polynomial.clone();
        }

        if self.randomness.len() == self.polynomial.number_of_vars() - 1 {
            self.randomness.push(Fp::sample());
            let latest_randomness = self.randomness.get(self.randomness.len() - 1).unwrap();
            assert_eq!(self.polynomial.eval(&*self.randomness), rec_polynomial.clone().unwrap().eval(&*vec![*latest_randomness]));
            self.interaction_completed = true;

            return;
        }

        self.step += 1;

        self.randomness.push(Fp::sample())
    }

    pub fn is_verifier_accept(&self) -> bool {
        self.interaction_completed
    }
}

#[cfg(test)]
mod tests {
    use crate::fp::Fp;
    use crate::mpolynomial::MPolynomial;
    use crate::SumcheckProtocol;

    #[test]
    fn test_sumcheck_protocol() {
        let coefficients = vec![Fp(2), Fp(3), Fp(5), Fp(7)];
        let powers = vec![1, 1, 1];

        let polynomial = MPolynomial::from(coefficients, powers);

        let mut protocol = SumcheckProtocol::new(polynomial);

        let step = protocol.prove();
        protocol.verify(step);

        let step = protocol.prove();
        protocol.verify(step);

        let step = protocol.prove();
        protocol.verify(step);

        assert!(protocol.is_verifier_accept())
    }
}