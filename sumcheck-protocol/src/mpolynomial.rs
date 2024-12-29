use crate::fp::Fp;
use crate::upolynomial::UPolynomial;

pub struct MPolynomial {
    pub coefficients: Vec<Fp>,
    pub powers: Vec<u32>,
}

impl MPolynomial {
    pub fn from(coefficients: Vec<Fp>, powers: Vec<u32>) -> Self {
        Self {
            coefficients,
            powers,
        }
    }

    pub fn eval(&self, x: &[Fp]) -> Fp {
        assert_eq!(x.len(), self.powers.len());

        let mut result = Fp(0);

        for i in 0..self.coefficients.len() {
            let mut term = self.coefficients[i];

            if i < self.powers.len() {
                result = result + (term * x[i].pow(self.powers[i]))
            } else {
                result = result + term;
            }
        }

        result
    }

    pub fn degree_ind(&self) -> usize {
        let mut result: usize = 0;

        for i in 0..self.powers.len() {
            if self.powers[i] > result as u32 {
                result = self.powers[i] as usize;
            }
        }

        result
    }

    pub fn sum_over_hyper_cube(&self, provided_x: Option<Vec<Fp>>) -> Fp {
        let vars_len = self.powers.len();
        let mut result = Fp(0);

        let mut initial_x = vec![];
        if provided_x.is_some() {
            initial_x = provided_x.unwrap();
        }

        for i in 0..2u64.pow(vars_len as u32 - initial_x.len() as u32) {
            let mut x = initial_x.clone();

            for j in 0..(vars_len - initial_x.len()) {
                x.push(Fp((i >> j) & 1));
            }

            result = result + self.eval(&x);
        }

        result
    }

    pub fn fix_var_over_hyper_cube(
        &self,
        provided_x: Option<&Vec<Fp>>,
    ) -> UPolynomial {
        let mut initial_x = vec![];
        if let Some(setup_x) = provided_x {
            initial_x = setup_x.clone();
        }

        let mut first_batch = initial_x.clone();
        let mut second_batch = initial_x.clone();

        first_batch.push(Fp(0));
        second_batch.push(Fp(1));

        let first_result = self.sum_over_hyper_cube(Some(first_batch));
        let second_result = self.sum_over_hyper_cube(Some(second_batch));

        UPolynomial::interpolate(vec![(Fp(0), first_result), (Fp(1), second_result)])
    }

    pub fn number_of_vars(&self) -> usize {
        self.powers.len()
    }
}

#[cfg(test)]
mod tests {
    use crate::fp::Fp;
    use crate::mpolynomial::MPolynomial;

    #[test]
    fn test_eval() {
        let coefficients = vec![Fp(3), Fp(5)];
        let powers = vec![2];
        let polynomial = MPolynomial::from(coefficients, powers);

        let x = vec![Fp(3)];
        assert_eq!(polynomial.eval(&x), Fp(3 * 3u64.pow(2) + 5))
    }

    #[test]
    fn test_sum_over_hyper_cube() {
        let coefficients = vec![Fp(2), Fp(3), Fp(5), Fp(7)];
        let powers = vec![1, 1, 1];

        let polynomial = MPolynomial::from(coefficients, powers);

        assert_eq!(polynomial.sum_over_hyper_cube(None), Fp(96))
    }

    #[test]
    fn test_sum_over_hyper_cube_with_args() {
        let coefficients = vec![Fp(2), Fp(3), Fp(5), Fp(7)];
        let powers = vec![1, 1, 1];

        let polynomial = MPolynomial::from(coefficients, powers);

        assert_eq!(polynomial.sum_over_hyper_cube(Some(vec![Fp(2)])), Fp(60))
    }
}
