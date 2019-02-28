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

#![cfg_attr(feature = "nightly", feature(specialization))]
pub extern crate cgmath;
extern crate num_traits;

#[macro_use]
#[cfg(feature="serde")]
extern crate serde;

macro_rules! d {
    ($($t:tt)*) => {
        <Self::D as Dimensionality<Self::Scalar>>::$($t)*
    };
}

#[cfg(feature="nightly")]
macro_rules! default_if_nightly {
    ($($t:tt)*) => {
        default $($t)*
    };
}

#[cfg(not(feature="nightly"))]
macro_rules! default_if_nightly {
    ($($t:tt)*) => {
        $($t)*
    };
}

pub mod ellipse;
pub mod line;
pub mod polar;
pub mod rect;

use cgmath::*;

pub trait Lerp: BaseScalarGeom {
    fn lerp<T: LerpFactor>(self, other: Self, t: T) -> Self {
        println!();
        Self::from(
            T::from(other).unwrap() * t +
            T::from(self).unwrap() * (T::one() - t)
        ).unwrap()
    }
}
impl<S: BaseScalarGeom> Lerp for S {}

pub trait LerpFactor: BaseScalarGeom + BaseFloat {}
impl<F: BaseScalarGeom + BaseFloat> LerpFactor for F {}

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
pub trait BasePointGeom:
          EuclideanSpace<Diff=<<Self as BasePointGeom>::D as Dimensionality<<Self as EuclideanSpace>::Scalar>>::Vector>
        + ElementWise<<Self as EuclideanSpace>::Scalar>
        + MulDiv<<Self as EuclideanSpace>::Scalar>
    where <Self::D as Dimensionality<<Self as EuclideanSpace>::Scalar>>::Vector: BaseVectorGeom<D=Self::D>,
          <Self::D as Dimensionality<<Self as EuclideanSpace>::Scalar>>::Vector: VectorSpace<Scalar=<Self as EuclideanSpace>::Scalar>,
          <Self as EuclideanSpace>::Scalar: BaseScalarGeom
{
    type D: Dimensionality<<Self as EuclideanSpace>::Scalar>;
}
pub trait BaseVectorGeom:
          VectorSpace
        + Array<Element=<Self as VectorSpace>::Scalar>
        + MulDiv
        + MulDiv<<Self as VectorSpace>::Scalar>
    where <Self as VectorSpace>::Scalar: BaseScalarGeom
{
    type D: Dimensionality<<Self as VectorSpace>::Scalar>;
}
impl<T: BaseNum + MulDiv + Bounded + AbsDistance> BaseScalarGeom for T {}
impl<S: BaseScalarGeom> BasePointGeom for Point1<S> { type D = D1; }
impl<S: BaseScalarGeom> BasePointGeom for Point2<S> { type D = D2; }
impl<S: BaseScalarGeom> BasePointGeom for Point3<S> { type D = D3; }
impl<S: BaseScalarGeom> BaseVectorGeom for Vector1<S> { type D = D1; }
impl<S: BaseScalarGeom> BaseVectorGeom for Vector2<S> { type D = D2; }
impl<S: BaseScalarGeom> BaseVectorGeom for Vector3<S> { type D = D3; }

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

pub trait Dimensionality<S: BaseScalarGeom>: 'static + Sized {
    type Point: BasePointGeom<D=Self> + EuclideanSpace<Scalar=S>;
    type Vector: BaseVectorGeom<D=Self> + VectorSpace<Scalar=S>;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub enum D1 {}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub enum D2 {}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub enum D3 {}
impl<S: BaseScalarGeom> Dimensionality<S> for D1 {
    type Point = Point1<S>;
    type Vector = Vector1<S>;
}
impl<S: BaseScalarGeom> Dimensionality<S> for D2 {
    type Point = Point2<S>;
    type Vector = Vector2<S>;
}
impl<S: BaseScalarGeom> Dimensionality<S> for D3 {
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
            default_if_nightly!{
                fn mul_div(mut self, mul: $array<S>, div: $array<S>) -> $array<S> {
                    for i in 0..$array::<S>::len() {
                        self[i] = self[i].mul_div(mul[i], div[i]);
                    }
                    self
                }
            }
        }

        impl<S: BaseNum + MulDiv> MulDiv<S> for $array<S> {
            #[inline(always)]
            default_if_nightly!{
                fn mul_div(mut self, mul: S, div: S) -> $array<S> {
                    for i in 0..$array::<S>::len() {
                        self[i] = self[i].mul_div(mul, div);
                    }
                    self
                }
            }
        }
    )*)
}

#[cfg(feature="nightly")]
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

        #[cfg(feature="nightly")]
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

#[cfg(feature="nightly")]
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

        #[cfg(feature = "nightly")]
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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn lerp() {
        assert_eq!(0, u32::lerp(0, 128, 0.0f32));
        assert_eq!(64, u32::lerp(0, 128, 0.5f32));
        assert_eq!(128, u32::lerp(0, 128, 1f32));

        assert_eq!(0, i32::lerp(i32::min_value(), i32::max_value(), 0.5f32));
        assert_eq!(u32::max_value() / 2, u32::lerp(0, u32::max_value(), 0.5f32));
    }
}
