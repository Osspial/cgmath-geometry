#![feature(specialization)]
extern crate cgmath;
extern crate num_traits;

pub mod rect;
pub mod line;

pub use self::rect::*;
pub use self::line::*;

use cgmath::*;
use num_traits::Bounded;

pub trait MulDiv<Rhs = Self> {
    fn mul_div(self, mul: Rhs, div: Rhs) -> Self;
}

pub trait BaseNumGeom: BaseNum + MulDiv + Bounded {}
impl<T: BaseNum + MulDiv + Bounded> BaseNumGeom for T {}

macro_rules! impl_mul_div_array_default {
    ($($array:ident),*) => ($(
        impl<S: BaseNum + MulDiv> MulDiv for $array<S> {
            #[inline(always)]
            default fn mul_div(mut self, mul: $array<S>, div: $array<S>) -> $array<S> {
                for i in 0..$array::<S>::len() {
                    self[i] = self[i].mul_div(mul[i], div[i]);
                }
                self
            }
        }

        impl<S: BaseNum + MulDiv> MulDiv<S> for $array<S> {
            #[inline(always)]
            default fn mul_div(mut self, mul: S, div: S) -> $array<S> {
                for i in 0..$array::<S>::len() {
                    self[i] = self[i].mul_div(mul, div);
                }
                self
            }
        }
    )*)
}

macro_rules! impl_mul_div_array_int {
    ($sm:ident => $bg:ident; $($array:ident),*) => {$(
        impl MulDiv for $array<$sm> {
            #[inline(always)]
            fn mul_div(self, mul: $array<$sm>, div: $array<$sm>) -> $array<$sm> {
                self.cast::<$bg>().mul_element_wise(mul.cast()).div_element_wise(div.cast()).cast()
            }
        }

        impl MulDiv<$sm> for $array<$sm> {
            #[inline(always)]
            fn mul_div(self, mul: $sm, div: $sm) -> $array<$sm> {
                self.cast::<$bg>().mul_element_wise(mul as $bg).div_element_wise(div as $bg).cast()
            }
        }
    )*};
}

macro_rules! impl_mul_div_int {
    ($(impl $sm:ident => $bg:ident;)*) => {$(
        impl MulDiv for $sm {
            #[inline(always)]
            fn mul_div(self, mul: $sm, div: $sm) -> $sm {
                let res = self as $bg * mul as $bg / div as $bg;
                debug_assert!(res <= $sm::max_value() as $bg);
                res as $sm
            }
        }

        impl_mul_div_array_int!($sm => $bg; Point1, Point2, Point3, Vector1, Vector2, Vector3, Vector4);
    )*}
}

macro_rules! impl_mul_div_array_float {
    ($float:ident; $($array:ident),*) => {$(
        impl MulDiv for $array<$float> {
            #[inline(always)]
            fn mul_div(self, mul: $array<$float>, div: $array<$float>) -> $array<$float> {
                self.mul_element_wise(mul).div_element_wise(div)
            }
        }

        impl MulDiv<$float> for $array<$float> {
            #[inline(always)]
            fn mul_div(self, mul: $float, div: $float) -> $array<$float> {
                self.mul_element_wise(mul).div_element_wise(div)
            }
        }
    )*};
}

macro_rules! impl_mul_div_float {
    ($(impl $float:ident;)*) => ($(
        impl MulDiv for $float {
            #[inline(always)]
            fn mul_div(self, mul: $float, div: $float) -> $float {
                self * mul / div
            }
        }

        impl_mul_div_array_float!($float; Point1);//, Point2, Point3, Vector1, Vector2, Vector3, Vector4);
    )*)
}

impl_mul_div_int! {
    impl u8 => u16;
    impl u16 => u32;
    impl u32 => u64;
    impl i8 => i16;
    impl i16 => i32;
    impl i32 => i64;
}

impl_mul_div_array_default!(Point1, Point2, Point3, Vector1, Vector2, Vector3, Vector4);

impl_mul_div_float! {
    impl f32;
    impl f64;
}
