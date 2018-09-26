// Copyright 2018 Osspial
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#![feature(specialization)]
pub extern crate cgmath;
extern crate num_traits;

#[macro_use]
#[cfg(feature="serde")]
extern crate serde;

macro_rules! d {
    ($($t:tt)*) => {
        <Self::D as Dimensionality>::$($t)*
    };
}

// pub mod ellipse;
pub mod line;
pub mod polar;
pub mod rect;

// pub use self::ellipse::*;
pub use self::rect::*;
pub use self::line::*;

use cgmath::*;
use std::ops::{Add, Sub, Mul, Div, Rem};
use std::iter::Sum;

pub trait MulDiv<Rhs = Self> {
    fn mul_div(self, mul: Rhs, div: Rhs) -> Self;
}

pub trait AbsDistance {
    type Abs: BaseScalarGeom;

    fn abs_distance(self, rhs: Self) -> Self::Abs;
    fn to_abs(self) -> Self::Abs;
    fn add_abs(self, rhs: Self::Abs) -> Self;
    fn sub_abs(self, rhs: Self::Abs) -> Self;
}

pub trait BaseScalarGeom: BaseNum + MulDiv + Bounded + AbsDistance {}
impl<T: BaseNum + MulDiv + Bounded + AbsDistance> BaseScalarGeom for T {}

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

pub trait Dimensionality {
    type Scalar: BaseScalarGeom;
    type Point:
          EuclideanSpace<Scalar=Self::Scalar, Diff=Self::Vector>
        + ElementWise<Self::Scalar>
        + MulDiv<Self::Scalar>;
    type Vector:
          VectorSpace<Scalar=Self::Scalar>
        + Array<Element=Self::Scalar>
        + MulDiv
        + MulDiv<Self::Scalar>;
}

pub struct D1<S: BaseScalarGeom>(S);
pub struct D2<S: BaseScalarGeom>(S);
pub struct D3<S: BaseScalarGeom>(S);
impl<S: BaseScalarGeom> Dimensionality for D1<S> {
    type Scalar = S;
    type Point = Point1<S>;
    type Vector = Vector1<S>;
}
impl<S: BaseScalarGeom> Dimensionality for D2<S> {
    type Scalar = S;
    type Point = Point2<S>;
    type Vector = Vector2<S>;
}
impl<S: BaseScalarGeom> Dimensionality for D3<S> {
    type Scalar = S;
    type Point = Point3<S>;
    type Vector = Vector3<S>;
}

impl<I> Into<Option<I>> for Intersection<I> {
    #[inline]
    fn into(self) -> Option<I> {
        match self {
            Intersection::Some(i) => Some(i),
            _ => None
        }
    }
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

macro_rules! impl_abs_dist_int {
    ($(impl $i:ty => $u:ty;)+) => {$(
        impl AbsDistance for $i {
            type Abs = $u;

            #[inline]
            fn abs_distance(self, rhs: $i) -> $u {
                let (min, max) = match self < rhs {
                    true => (self, rhs),
                    false => (rhs, self)
                };
                max.wrapping_sub(min) as $u
            }

            #[inline]
            #[allow(unused_comparisons)]
            fn to_abs(self) -> $u {
                if self < 0 {
                    (0 as $i).wrapping_sub(self) as $u
                } else {
                    self as $u
                }
            }

            #[inline]
            fn add_abs(self, rhs: $u) -> $i {
                self.wrapping_add(rhs as $i)
            }

            #[inline]
            fn sub_abs(self, rhs: $u) -> $i {
                self.wrapping_sub(rhs as $i)
            }
        }
    )+};
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

        impl AbsDistance for $float {
            type Abs = $float;
            #[inline]
            fn abs_distance(self, rhs: $float) -> $float {
                let (min, max) = match self < rhs {
                    true => (self, rhs),
                    false => (rhs, self)
                };
                max - min
            }

            #[inline]
            fn to_abs(self) -> $float {
                self.abs()
            }

            #[inline]
            fn add_abs(self, rhs: $float) -> $float {
                self + rhs
            }

            #[inline]
            fn sub_abs(self, rhs: $float) -> $float {
                self - rhs
            }
        }

        impl_mul_div_array_float!($float; Point1, Point2, Point3, Vector1, Vector2, Vector3, Vector4);
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

impl_abs_dist_int! {
    impl i8 => u8;
    impl u8 => u8;
    impl i16 => u16;
    impl u16 => u16;
    impl i32 => u32;
    impl u32 => u32;
    // impl i64 => u64;
    // impl u64 => u64;
}

impl_mul_div_float! {
    impl f32;
    impl f64;
}

impl_mul_div_array_default!(Point1, Point2, Point3, Vector1, Vector2, Vector3, Vector4);
