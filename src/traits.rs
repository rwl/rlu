use std::{fmt, ops};

use num_complex::Complex;

pub trait Int: num_traits::PrimInt + fmt::Display + fmt::Debug {
    fn from_usize(i: usize) -> Self {
        match Self::from(i) {
            Some(j) => j,
            None => panic!("must be able to create Int from {}", i),
        }
    }

    fn to_index(&self) -> usize {
        match self.to_usize() {
            Some(j) => j,
            None => panic!("must be able to convert Int to usize: {}", self),
        }
    }
}

impl Int for usize {}
impl Int for u8 {}
impl Int for u16 {}
impl Int for u32 {}
impl Int for u64 {}
impl Int for u128 {}

impl Int for isize {}
impl Int for i8 {}
impl Int for i16 {}
impl Int for i32 {}
impl Int for i64 {}
impl Int for i128 {}

pub trait Scalar:
    Copy
    // + PartialOrd
    + num_traits::Zero
    + ops::Mul<Output = Self>
    + ops::Div<Output = Self>
    + ops::SubAssign
    + ops::DivAssign
    + Norm<Self::Norm>
    + fmt::Display
    + fmt::Debug
{
    type Norm: PartialOrd + fmt::Display;

    #[cfg(feature = "debug")]
    fn pretty_string(&self, _config: pretty_dtoa::FmtFloatConfig) -> String {
        format!("{}", self)
    }
}

impl Scalar for f64 {
    type Norm = f64;

    #[cfg(feature = "debug")]
    fn pretty_string(&self, config: pretty_dtoa::FmtFloatConfig) -> String {
        pretty_dtoa::dtoa(*self, config)
    }
}

impl Scalar for Complex<f64> {
    type Norm = f64;

    #[cfg(feature = "debug")]
    fn pretty_string(&self, config: pretty_dtoa::FmtFloatConfig) -> String {
        format!(
            "{}{}j{}",
            pretty_dtoa::dtoa(self.re, config),
            if self.im.signum() < 0.0 { "-" } else { "+" },
            pretty_dtoa::dtoa(self.im, config)
        )
        .to_string()
    }
}

impl Scalar for usize {
    type Norm = Self;
}
// TODO: ints

pub trait Norm<F> {
    fn norm(&self) -> F;
}

impl Norm<f64> for f64 {
    fn norm(&self) -> f64 {
        f64::abs(*self)
    }
}

impl Norm<f32> for f32 {
    fn norm(&self) -> f32 {
        f32::abs(*self)
    }
}

impl Norm<f64> for Complex<f64> {
    fn norm(&self) -> f64 {
        num_complex::Complex::norm(*self)
    }
}

impl Norm<usize> for usize {
    fn norm(&self) -> usize {
        *self
    }
}
// TODO: ints
