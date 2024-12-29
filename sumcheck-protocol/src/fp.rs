use rand::random;
use std::cmp::PartialEq;
use std::ops::{Add, Div, Mul, Neg, Rem, Shr, ShrAssign, Sub, SubAssign};

pub trait FiniteField {
    /// Goldilocks prime under u64 (dec: 18446744069414584321)
    const MODULO: u64 = (2u128.pow(64) - 2u128.pow(32) + 1) as u64;
}

/// Could be optimized using Montgomery multiplication
/// See section 2.3.2 of Tolga Acar's thesis
/// https://www.microsoft.com/en-us/research/wp-content/uploads/1998/06/97Acar.pdf
#[derive(Debug, Clone, Copy)]
pub struct Fp(pub u64);

impl Fp {
    pub fn neg(&self) -> Self {
        Fp(Self::MODULO - self.0)
    }

    pub fn pow(&self, p0: u32) -> Self {
        let mut result = Fp(1);
        let mut base = *self;
        let mut p = p0;

        while p > 0 {
            if p & 1 == 1 {
                result = result * base;
            }

            base = base * base;
            p >>= 1;
        }

        result
    }
}

impl Fp {
    pub const MAX: Self = Fp(Fp::MODULO - 1);

    pub fn zero() -> Self {
        Fp(0)
    }

    pub fn one() -> Self {
        Fp(1)
    }

    pub fn sample() -> Self {
        Fp::from(random::<u64>())
    }

    pub fn from<T: Into<u128>>(value: T) -> Self {
        Fp((value.into() % Self::MODULO as u128) as u64)
    }

    /// Algorithm 16 in "Efficient Software-Implementation of Finite Fields with Applications to Cryptography"
    pub fn inverse(&self) -> Self {
        let modulo = Self::MODULO;

        if self.0 == 0 {
            return Fp(0);
        }

        let mut v: u64 = self.0;
        let mut u: u64 = Self::MODULO;
        let mut s: u64 = 1;
        let mut r: u64 = 0;

        let mut carry: bool;
        let mut borrow: bool;

        while u != 1 && v != 1 {
            while v & 1 == 0 {
                v >>= 1;
                if s & 1 == 0 {
                    s >>= 1;
                } else {
                    (s, carry) = u64::overflowing_add(s, modulo);
                    s >>= 1;
                    if carry {
                        s |= 1 << 63;
                    }
                }
            }

            while u & 1 == 0 {
                u >>= 1;
                if r & 1 == 0 {
                    r >>= 1;
                } else {
                    (r, carry) = u64::overflowing_add(r, modulo);
                    r >>= 1;
                    if carry {
                        r |= 1 << 63;
                    }
                }
            }

            if v >= u {
                v -= u;
                (s, borrow) = u64::overflowing_sub(s, r);
                if borrow {
                    s = u64::overflowing_add(s, modulo).0;
                }
            } else {
                u -= v;
                (r, borrow) = u64::overflowing_sub(r, s);
                if borrow {
                    r = u64::overflowing_add(r, modulo).0;
                }
            }
        }

        if u == 1 {
            Fp::from(r)
        } else {
            Fp::from(s)
        }
    }

    pub fn multi_inv(a: &Vec<Fp>) -> Vec<Fp> {
        let mut partials = vec![Fp::zero(); a.len() + 1];

        partials[0] = Fp::one();
        for i in 0..a.len() {
            partials[i + 1] = partials[i] * a[i];
        }

        let mut inv = Fp::inverse(&partials[a.len()]);
        let mut outputs = vec![Fp::zero(); a.len()];

        for i in (0..a.len()).rev() {
            outputs[i] = partials[i] * inv;

            if a[i] == Fp::zero() {
                outputs[i] = Fp::one();
            }

            inv = inv * a[i];
        }

        outputs
    }
}

impl FiniteField for Fp {}
impl FiniteField for &Fp {}

impl Add for Fp {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Fp(((self.0 as u128 + rhs.0 as u128) % Self::MODULO as u128) as u64)
    }
}

impl Sub for Fp {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Fp(((self.0 as u128 + Self::MODULO as u128 - rhs.0 as u128) % Self::MODULO as u128) as u64)
    }
}

impl SubAssign for Fp {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Mul for Fp {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Fp(((self.0 as u128 * rhs.0 as u128) % Self::MODULO as u128) as u64)
    }
}

impl Div for Fp {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        // Sanitize the input
        let a = Fp::from(self.0);
        let b = Fp::from(rhs.0);

        a * b.inverse()
    }
}

impl Rem for Fp {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        Fp(((self.0 as u128 % rhs.0 as u128) % Self::MODULO as u128) as u64)
    }
}

impl Shr<u64> for Fp {
    type Output = Self;

    fn shr(self, rhs: u64) -> Self::Output {
        Fp(self.0 >> rhs)
    }
}

impl ShrAssign<u64> for Fp {
    fn shr_assign(&mut self, rhs: u64) {
        self.0 >>= rhs;
    }
}

impl Neg for Fp {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Fp(Self::MODULO - self.0)
    }
}

impl PartialEq for Fp {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }

    fn ne(&self, other: &Self) -> bool {
        self.0 != other.0
    }
}

impl PartialEq<Fp> for &Fp {
    fn eq(&self, other: &Fp) -> bool {
        self.0 == other.0
    }
}

impl PartialOrd for Fp {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

#[cfg(test)]
mod tests {
    use crate::fp::{FiniteField, Fp};

    #[test]
    fn test_add() {
        let a = Fp(120);
        let b = Fp(80);

        assert_eq!(a + b, Fp(200))
    }

    #[test]
    fn test_add_max() {
        let a = Fp(u64::MAX);
        let b = Fp(u64::MAX);

        assert_eq!(a + b, Fp(8589934588))
    }

    #[test]
    fn test_sub() {
        let a = Fp(120);
        let b = Fp(80);

        assert_eq!(b - a, Fp(Fp::MODULO - 40))
    }

    #[test]
    fn test_sub_max() {
        let a = Fp(u64::MAX);
        let b = Fp(u64::MAX);

        assert_eq!(a - b, Fp(0))
    }

    #[test]
    fn test_mul() {
        let a = Fp(120);
        let b = Fp(80);

        assert_eq!(a * b, Fp(120 * 80))
    }

    #[test]
    fn test_mul_max() {
        let a = Fp(u64::MAX);
        let b = Fp(u64::MAX);

        assert_eq!(a * b, Fp(18446744056529682436))
    }

    #[test]
    fn test_div() {
        let a = Fp(120);
        let b = Fp(4);

        assert_eq!(a / b, Fp(120 / 4))
    }

    #[test]
    fn test_div_max() {
        let a = Fp(u64::MAX);
        let b = Fp(u64::MAX);

        assert_eq!(a / b, Fp(1))
    }

    #[test]
    fn test_pow() {
        let a = Fp(2);
        let b = 10;

        assert_eq!(a.pow(b), Fp(1024))
    }
}
