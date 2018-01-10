#![feature(specialization)]
extern crate cgmath;
extern crate num_traits;

macro_rules! P {
    ($($t:tt)*) => (<P as EuclideanSpace>$($t)*);
}

pub mod ellipse;
pub mod line;
pub mod polar;
pub mod rect;

pub use self::ellipse::*;
pub use self::rect::*;
pub use self::line::*;

use cgmath::*;
use num_traits::Bounded;

pub trait MulDiv<Rhs = Self> {
    fn mul_div(self, mul: Rhs, div: Rhs) -> Self;
}

pub trait BaseScalarGeom: BaseNum + MulDiv + Bounded {}
impl<T: BaseNum + MulDiv + Bounded> BaseScalarGeom for T {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Intersection<I> {
    Some(I),
    Eq,
    None
}

pub trait Intersect<RHS=Self> {
    type Intersection;
    fn intersect(self, rhs: RHS) -> Intersection<Self::Intersection>;
}

fn cmp_min<S: BaseNum>(l: S, r: S) -> S {
    if l < r {
        l
    } else {
        r
    }
}

fn cmp_max<S: BaseNum>(l: S, r: S) -> S {
    if l > r {
        l
    } else {
        r
    }
}

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
    ($sm:ident => $bg:ident; $($array:ident{ $($field:ident),+ }),*) => {$(
        impl MulDiv for $array<$sm> {
            #[inline(always)]
            fn mul_div(self, mul: $array<$sm>, div: $array<$sm>) -> $array<$sm> {
                fn cast_no_check(sm: $array<$sm>) -> $array<$bg> {
                    $array{ $($field: sm.$field as $bg),+ }
                }
                let res = cast_no_check(self).mul_element_wise(cast_no_check(mul)).div_element_wise(cast_no_check(div));
                $(
                    debug_assert!(res.$field <= $sm::max_value() as $bg);
                )+
                $array{ $($field: res.$field as $sm),+ }
            }
        }

        impl MulDiv<$sm> for $array<$sm> {
            #[inline(always)]
            fn mul_div(self, mul: $sm, div: $sm) -> $array<$sm> {
                fn cast_no_check(sm: $array<$sm>) -> $array<$bg> {
                    $array{ $($field: sm.$field as $bg),+ }
                }
                let res = cast_no_check(self).mul_element_wise(mul as $bg).div_element_wise(div as $bg);
                $(
                    debug_assert!(res.$field <= $sm::max_value() as $bg);
                )+
                $array{ $($field: res.$field as $sm),+ }
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

        impl_mul_div_array_int!($sm => $bg; Point1{x}, Point2{x, y}, Point3{x, y, z}, Vector1{x}, Vector2{x, y}, Vector3{x, y, z}, Vector4{x, y, z, w});
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
