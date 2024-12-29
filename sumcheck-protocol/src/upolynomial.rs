use crate::fp::Fp;
use std::ops::Div;

#[derive(Debug, Clone)]
pub struct UPolynomial {
    pub coefficients: Vec<Fp>,
}

impl UPolynomial {
    pub fn from(coefficients: Vec<Fp>) -> Self {
        Self { coefficients }
    }

    pub fn zero_at_given_x(xs: &Vec<Fp>) -> Self {
        let mut root = vec![Fp::one()];

        for x in xs {
            root.insert(0, Fp::zero());

            for j in 0..root.len() - 1 {
                root[j] = root[j] - (root[j + 1] * *x)
            }
        }

        UPolynomial::from(root)
    }

    pub fn interpolate(points: Vec<(Fp, Fp)>) -> Self {
        let (xs, ys) = &points
            .iter()
            .map(|(x, y)| (x.clone(), y.clone()))
            .unzip::<Fp, Fp, Vec<Fp>, Vec<Fp>>();

        let mut root = UPolynomial::zero_at_given_x(xs);

        if root.coefficients.len() != points.len() + 1 {
            return UPolynomial::from(vec![]);
        }

        let mut numerators = vec![];
        for x in xs {
            numerators.push(root.clone() / UPolynomial::from(vec![Fp::neg(x), Fp::one()]))
        }

        let mut denominator = vec![];
        for i in 0..xs.len() {
            denominator.push(numerators[i].eval(&[xs[i]]))
        }

        let inv_denominators = Fp::multi_inv(&denominator);

        let mut b = vec![];
        for i in 0..xs.len() {
            let y_slice = ys[i] * inv_denominators[i];

            for j in 0..ys.len() {
                if numerators[i].coefficients[j] != Fp::zero() && ys[i] != Fp::zero() {
                    let res_mul = numerators[i].coefficients[j] * y_slice;

                    if b.len() <= j {
                        b.push(Fp::zero());
                    }

                    b[j] = b[j] + res_mul;
                }
            }
        }

        UPolynomial::from(b)
    }

    pub fn eval(&self, x: &[Fp]) -> Fp {
        let coefficients_len = self.coefficients.len();

        assert_eq!(x.len(), coefficients_len - 1);

        let mut result = Fp(0);

        for i in 0..coefficients_len {
            let mut term = self.coefficients[i];

            if i < x.len() {
                result = result + (term * x[i].pow(coefficients_len as u32 - 1 - i as u32))
            } else {
                result = result + term;
            }
        }

        result
    }

    pub fn degree(&self) -> usize {
        self.coefficients.len() - 1
    }
}

impl Div for UPolynomial {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let mut a_copy = self.coefficients.clone();

        if a_copy.len() < rhs.coefficients.len() {
            return UPolynomial::from(vec![]);
        }

        let mut output = vec![];

        let mut a_pos = a_copy.len() - 1;
        let mut b_pos = rhs.coefficients.len() - 1;
        let mut diff = a_pos - b_pos;

        loop {
            let quot = a_copy[a_pos] / rhs.coefficients[b_pos];
            output.push(quot);

            for i in (0..=b_pos).rev() {
                a_copy[diff + i] = a_copy[diff + i] - (rhs.coefficients[i] * quot);
            }

            a_pos -= 1;

            if diff == 0 {
                break;
            }

            diff -= 1;
        }

        UPolynomial::from(output)
    }
}

#[cfg(test)]
mod tests {
    use crate::fp::Fp;
    use crate::upolynomial::UPolynomial;

    #[test]
    fn test_eval() {
        let coefficients = vec![Fp(8), Fp(3), Fp(5)];
        let polynomial = UPolynomial::from(coefficients);

        let x = vec![Fp(3), Fp(5)];
        assert_eq!(
            polynomial.eval(&x),
            Fp(8 * 3u64.pow(2) + 3 * 5u64.pow(1) + 5)
        )
    }

    #[test]
    fn test_interpolate() {
        let points = vec![(Fp(0), Fp(44)), (Fp(1), Fp(52))];
        let polynomial = UPolynomial::interpolate(points);

        assert_eq!(polynomial.coefficients, vec![Fp(8), Fp(44)])
    }
}
