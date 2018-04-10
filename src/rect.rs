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

use {MulDiv, BaseScalarGeom, AbsDistance};
use cgmath::*;

use line::{Linear, Segment};

use std::ops::{Add, Sub};
use num_traits::{Bounded, NumCast, ToPrimitive};

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct DimsBox<P: EuclideanSpace> {
    pub dims: P::Diff
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct OffsetBox<P: EuclideanSpace> {
    pub origin: P,
    pub dims: P::Diff
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct BoundBox<P: EuclideanSpace> {
    pub min: P,
    pub max: P
}

pub trait GeoBox {
    type Scalar: BaseScalarGeom;
    type Point: EuclideanSpace<Scalar=Self::Scalar, Diff=Self::Vector> + ElementWise<Self::Scalar> + MulDiv<Self::Scalar>;
    type Vector: VectorSpace<Scalar=Self::Scalar> + Array<Element=Self::Scalar> + MulDiv + MulDiv<Self::Scalar> + ElementWise;

    fn from_bounds(min: Self::Point, max: Self::Point) -> Self;

    #[inline]
    fn min(&self) -> Self::Point {
        self.max() - self.dims()
    }

    #[inline]
    fn max(&self) -> Self::Point {
        self.min() + self.dims()
    }

    #[inline]
    fn dims(&self) -> Self::Vector {
        self.max() - self.min()
    }

    #[inline]
    fn width(&self) -> Self::Scalar {
        match Self::Point::len() >= 1 {
            true => self.dims()[0],
            false => Self::Scalar::zero()
        }
    }
    #[inline]
    fn height(&self) -> Self::Scalar {
        match Self::Point::len() >= 2 {
            true => self.dims()[1],
            false => Self::Scalar::zero()
        }
    }
    #[inline]
    fn depth(&self) -> Self::Scalar {
        match Self::Point::len() >= 3 {
            true => self.dims()[2],
            false => Self::Scalar::zero()
        }
    }

    #[inline]
    fn center(&self) -> Self::Point {
        self.min() + (self.dims() / (Self::Scalar::one() + Self::Scalar::one()))
    }

    #[inline]
    fn contains(&self, point: Self::Point) -> bool {
        let min = self.min();
        let max = self.max();

        let mut contains = true;
        for i in 0..Self::Point::len() {
            contains = contains &&
                min[i] <= point[i] &&
                point[i] <= max[i];
        }

        contains
    }

    fn intersect_rect(&self, other: Self) -> Option<Self>
        where Self: Sized
    {
        let (s_min, s_max) = (self.min(), self.max());
        let (o_min, o_max) = (other.min(), other.max());

        let (mut min, mut max) = (s_min, s_max);

        for i in 0..Self::Point::len() {
            min[i] = ::cmp_max(s_min[i], o_min[i]);
            max[i] = ::cmp_min(s_max[i], o_max[i]);

            if max[i] < min[i] {
                return None;
            }
        }

        Some(Self::from_bounds(min, max))
    }

    fn intersect_line<L>(&self, line: L) -> (Option<Self::Point>, Option<Self::Point>)
        where L: Linear<Scalar=Self::Scalar, Point=Self::Point, Vector=Self::Vector>,
              Self::Scalar: Bounded
    {
        macro_rules! switch_fn {
            ($name:ident) => {
                fn $name<R, L>(rect: &R, line: L) -> (Option<R::Point>, Option<R::Point>)
                    where R: GeoBox<Scalar=S> + ?Sized,
                          R::Point: EuclideanSpace<Scalar=S, Diff=R::Vector> + ElementWise<S> + MulDiv<S>,
                          R::Vector: VectorSpace<Scalar=S> + Array<Element=S> + MulDiv + MulDiv<S>,
                          L: Linear<Scalar=R::Scalar, Point=R::Point, Vector=R::Vector>;
            };
            ($({$prefix:tt})* $name:ident($rect:ident, $line:ident), $body:block) => {
                $($prefix)* fn $name<R, L>($rect: &R, $line: L) -> (Option<R::Point>, Option<R::Point>)
                    where R: GeoBox<Scalar=S> + ?Sized,
                          R::Point: EuclideanSpace<Scalar=S, Diff=R::Vector> + ElementWise<S> + MulDiv<S>,
                          R::Vector: VectorSpace<Scalar=S> + Array<Element=S> + MulDiv + MulDiv<S>,
                          L: Linear<Scalar=R::Scalar, Point=R::Point, Vector=R::Vector>
                $body
            }
        }
        struct TS;
        trait TypeSwitch<S>
            where S: BaseNum + MulDiv + Bounded + AbsDistance
        {
            switch_fn!{intersect_ts}
        }

        impl<S> TypeSwitch<S> for TS
            where S: BaseFloat + BaseNum + MulDiv + Bounded + AbsDistance
        {
            switch_fn!{intersect_ts(rect, line), {
                let line_origin = line.origin();
                let dir = line.dir();
                let dir_recip = line.dir_recip();
                let (rect_min, rect_max) = (rect.min(), rect.max());

                let (mut t_min, mut t_max) = (S::neg_infinity(), S::infinity());

                for i in 0..L::Point::len() {
                    let t_enter = (rect_min[i] - line_origin[i]) * dir_recip[i];
                    let t_exit = (rect_max[i] - line_origin[i]) * dir_recip[i];
                    t_min = t_min.max(t_enter.min(t_exit));
                    t_max = t_max.min(t_enter.max(t_exit));
                }

                if t_max < t_min {
                    (None, None)
                } else {
                    let t_of_point = |point: L::Point| {
                        let mut t = S::zero();
                        for i in 0..L::Point::len() {
                            let t_axis = (point[i] - line_origin[i]) * dir_recip[i];
                            if t_axis.abs() > t.abs() {
                                t = t_axis;
                            }
                        }
                        t
                    };
                    let t_enter = match line.start() {
                        None => Some(t_min),
                        Some(start) => match t_min >= t_of_point(start) {
                            true => Some(t_min),
                            false => None
                        }
                    };
                    let t_exit = match line.end() {
                        None => Some(t_max),
                        Some(end) => match t_max <= t_of_point(end) {
                            true => Some(t_max),
                            false => None
                        }
                    };
                    (t_enter.map(|t| line_origin + dir * t), t_exit.map(|t| line_origin + dir * t))
                }
            }}
        }

        impl<S> TypeSwitch<S> for TS
            where S: BaseNum + MulDiv + Bounded + AbsDistance
        {
            switch_fn!{{default} intersect_ts(rect, line), {
                let zero = L::Scalar::zero();
                let min = rect.min();
                let max = rect.max();

                let Segment{ start, end } = line.clip_to_scalar_bounds();
                let (mut enter, mut exit) = (start, end);
                let (mut enter_valid, mut exit_valid) = (false, false);
                let dir = line.dir();

                for i in 0..R::Point::len() {
                    if enter[i] <= exit[i] {
                        if enter[i] <= min[i] && min[i] <= exit[i] && dir[i] != zero {
                            enter = enter + dir.mul_div(min[i] - enter[i], dir[i]);
                            enter_valid = true;
                        }

                        if enter[i] <= max[i] && max[i] < exit[i] && dir[i] != zero {
                            exit = exit - dir.mul_div(exit[i] - max[i], dir[i]);
                            exit_valid = true;
                        }
                    } else {
                        if exit[i] <= max[i] && max[i] <= enter[i] && dir[i] != zero {
                            enter = enter - dir.mul_div(enter[i] - max[i], dir[i]);
                            enter_valid = true;
                        }

                        if exit[i] < min[i] && min[i] <= enter[i] && dir[i] != zero {
                            exit = exit + dir.mul_div(min[i] - exit[i], dir[i]);
                            exit_valid = true;
                        }
                    };
                }

                for i in 0..R::Point::len() {
                    if enter[i] < min[i] || max[i] < enter[i] {
                        enter_valid = false;
                    }
                    if exit[i] < min[i] || max[i] < exit[i] {
                        exit_valid = false;
                    }
                }

                (
                    match enter_valid {
                        true => Some(enter),
                        false => None
                    },
                    match exit_valid {
                        true => Some(exit),
                        false => None
                    }
                )
            }}
        }

        TS::intersect_ts(self, line)
    }
}

macro_rules! inherent_impl_dims_offset {
    ($PointN:ident, $VectorN:ident; $new:ident; $($origin:ident, $dim:ident),+) => {
        impl<S: BaseScalarGeom> DimsBox<$PointN<S>> {
            #[inline]
            pub fn $new($($dim: S),+) -> DimsBox<$PointN<S>> {
                DimsBox {
                    dims: $VectorN::new($($dim),+)
                }
            }

            #[inline]
            pub fn cast<T>(&self) -> Option<DimsBox<$PointN<T>>>
                where T: NumCast + BaseScalarGeom,
                      S: ToPrimitive
            {
                Some(DimsBox {
                    dims: $VectorN::cast(&self.dims)?
                })
            }
        }

        impl<S: BaseScalarGeom> OffsetBox<$PointN<S>> {
            #[inline]
            pub fn $new($($origin: S,)+ $($dim: S),+) -> OffsetBox<$PointN<S>> {
                OffsetBox {
                    origin: $PointN::new($($origin),+),
                    dims: $VectorN::new($($dim),+)
                }
            }

            #[inline]
            pub fn cast<T>(&self) -> Option<OffsetBox<$PointN<T>>>
                where T: NumCast + BaseScalarGeom,
                      S: ToPrimitive
            {
                Some(OffsetBox {
                    origin: $PointN::cast(&self.origin)?,
                    dims: $VectorN::cast(&self.dims)?
                })
            }
        }
    }
}

macro_rules! inherent_impl_bounds {
    ($PointN:ident, $VectorN:ident; $new:ident; ($($min:ident),+), ($($max:ident),+)) => {
        impl<S: BaseScalarGeom> BoundBox<$PointN<S>> {
            #[inline]
            pub fn $new($($min: S),+, $($max: S),+) -> BoundBox<$PointN<S>> {
                BoundBox {
                    min: $PointN::new($($min),+),
                    max: $PointN::new($($max),+)
                }
            }

            #[inline]
            pub fn cast<T>(&self) -> Option<BoundBox<$PointN<T>>>
                where T: NumCast + BaseScalarGeom,
                      S: ToPrimitive
            {
                Some(BoundBox {
                    min: $PointN::cast(&self.min)?,
                    max: $PointN::cast(&self.max)?,
                })
            }
        }
    }
}

impl<P: EuclideanSpace> DimsBox<P> {
    #[inline]
    pub fn new(dims: P::Diff) -> DimsBox<P> {
        DimsBox{ dims }
    }
}
impl<P: EuclideanSpace> OffsetBox<P> {
    #[inline]
    pub fn new(origin: P, dims: P::Diff) -> OffsetBox<P> {
        OffsetBox{ origin, dims }
    }
}
impl<P: EuclideanSpace> BoundBox<P> {
    #[inline]
    pub fn new(min: P, max: P) -> BoundBox<P> {
        BoundBox{ min, max }
    }
}

inherent_impl_dims_offset!(Point1, Vector1; new1; origin_x, width);
inherent_impl_dims_offset!(Point2, Vector2; new2; origin_x, width, origin_y, height);
inherent_impl_dims_offset!(Point3, Vector3; new3; origin_x, width, origin_y, height, origin_z, depth);

inherent_impl_bounds!(Point1, Vector1; new1; (min_x), (max_x));
inherent_impl_bounds!(Point2, Vector2; new2; (min_x, min_y), (max_x, max_y));
inherent_impl_bounds!(Point3, Vector3; new3; (min_x, min_y, min_z), (max_x, max_y, max_z));

impl<P> GeoBox for DimsBox<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise
{
    type Scalar = P!(::Scalar);
    type Point = P;
    type Vector = P::Diff;

    #[inline]
    fn from_bounds(min: P, max: P) -> DimsBox<P> {
        DimsBox {
            dims: max - min
        }
    }

    #[inline]
    fn min(&self) -> P {P::from_value(P::Scalar::zero())}
    #[inline]
    fn dims(&self) -> P::Diff {self.dims}
}

impl<P> Bounded for DimsBox<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise + Bounded
{
    #[inline]
    fn min_value() -> DimsBox<P> {
        DimsBox {
            dims: P::Diff::min_value()
        }
    }

    #[inline]
    fn max_value() -> DimsBox<P> {
        DimsBox {
            dims: P::Diff::max_value()
        }
    }
}

impl<P> GeoBox for OffsetBox<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise
{
    type Scalar = P::Scalar;
    type Point = P;
    type Vector = P::Diff;

    #[inline]
    fn from_bounds(min: P, max: P) -> OffsetBox<P> {
        OffsetBox {
            origin: min,
            dims: max - min
        }
    }

    #[inline]
    fn min(&self) -> P {self.origin}
    #[inline]
    fn dims(&self) -> P::Diff {self.dims}
}

impl<P> GeoBox for BoundBox<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise
{
    type Scalar = P::Scalar;
    type Point = P;
    type Vector = P::Diff;

    #[inline]
    fn from_bounds(min: P, max: P) -> BoundBox<P> {
        BoundBox {
            min, max
        }
    }

    #[inline]
    fn min(&self) -> P {self.min}
    #[inline]
    fn max(&self) -> P {self.max}
}

impl<P> From<DimsBox<P>> for OffsetBox<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise
{
    #[inline]
    fn from(rect: DimsBox<P>) -> OffsetBox<P> {
        OffsetBox {
            origin: P::from_value(P::Scalar::zero()),
            dims: rect.dims
        }
    }
}

impl<P> From<DimsBox<P>> for BoundBox<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise
{
    #[inline]
    fn from(rect: DimsBox<P>) -> BoundBox<P> {
        BoundBox {
            min: P::from_value(P::Scalar::zero()),
            max: P::from_vec(rect.dims)
        }
    }
}

impl<P> From<OffsetBox<P>> for BoundBox<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise
{
    #[inline]
    fn from(rect: OffsetBox<P>) -> BoundBox<P> {
        BoundBox {
            min: rect.origin,
            max: rect.origin + rect.dims
        }
    }
}

impl<P> From<BoundBox<P>> for OffsetBox<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise
{
    #[inline]
    fn from(rect: BoundBox<P>) -> OffsetBox<P> {
        OffsetBox {
            origin: rect.min,
            dims: rect.max.to_vec() - rect.min.to_vec()
        }
    }
}

impl<P> Add<P::Diff> for OffsetBox<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise
{
    type Output = Self;
    #[inline]
    fn add(mut self, rhs: P::Diff) -> OffsetBox<P> {
        self.origin = self.origin + rhs;
        self
    }
}

impl<P> Sub<P::Diff> for OffsetBox<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise
{
    type Output = Self;
    #[inline]
    fn sub(mut self, rhs: P::Diff) -> OffsetBox<P> {
        self.origin = self.origin - rhs;
        self
    }
}

impl<P> Add<P::Diff> for BoundBox<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise
{
    type Output = Self;
    #[inline]
    fn add(mut self, rhs: P::Diff) -> BoundBox<P> {
        self.min = self.min + rhs;
        self.max = self.max + rhs;
        self
    }
}

impl<P> Sub<P::Diff> for BoundBox<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise
{
    type Output = Self;
    #[inline]
    fn sub(mut self, rhs: P::Diff) -> BoundBox<P> {
        self.min = self.min - rhs;
        self.max = self.max - rhs;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use line::Segment;

    #[test]
    fn line_touch_edge_exit() {
        let rect = BoundBox::new2(20, 20, 30, 40);
        let segment = Segment::new2(25, 25, 20, 25);
        assert_eq!((None, None), rect.intersect_line(segment));
    }

    #[test]
    fn line_touch_edge_outer_exit() {
        let rect = BoundBox::new2(20, 20, 30, 40);
        let segment = Segment::new2(25, 15, 25, 20);
        assert_eq!((Some(Point2::new(25, 20)), None), rect.intersect_line(segment));
    }

    #[test]
    fn line_exit_single_move_lt() {
        let rect = BoundBox::new2(20, 20, 30, 40);
        let segment = Segment::new2(20, 25, 19, 25);
        assert_eq!((None, Some(Point2::new(20, 25))), rect.intersect_line(segment));
    }

    #[test]
    fn line_exit_single_move_gt() {
        let rect = BoundBox::new2(20, 20, 30, 40);
        let segment = Segment::new2(30, 25, 31, 25);
        assert_eq!((None, Some(Point2::new(30, 25))), rect.intersect_line(segment));
    }

    #[test]
    fn line_touch_edge_enter() {
        let rect = BoundBox::new2(20, 20, 30, 40);
        let segment = Segment::new2(25, 0, 27, 20);
        assert_eq!((Some(Point2::new(27, 20)), None), rect.intersect_line(segment));
    }

    #[test]
    fn intersect_diagonal() {
        let rect = BoundBox::new2(20, 20, 40, 40);
        let segment = Segment::new2(0, 0, 50, 50);
        assert_eq!((Some(Point2::new(20, 20)), Some(Point2::new(40, 40))), rect.intersect_line(segment));
    }

    #[test]
    fn intersect_horizontal() {
        let rect = BoundBox::new2(20, 20, 30, 40);
        let segment = Segment::new2(0, 25, 50, 25);
        assert_eq!((Some(Point2::new(20, 25)), Some(Point2::new(30, 25))), rect.intersect_line(segment));
    }

    #[test]
    fn test_intersect() {
        assert_eq!(BoundBox::new2(20., 20., 40., 40.).intersect_rect(BoundBox::new2(0., 0., 10., 10.)), None);
        assert_eq!(BoundBox::new2(0., 0., 10., 10.).intersect_rect(BoundBox::new2(20., 20., 40., 40.)), None);
        assert_eq!(BoundBox::new2(10., 10., 20., 20.).intersect_rect(BoundBox::new2(5., 5., 15., 15.)), Some(BoundBox::new2(10., 10., 15., 15.)));
    }
}
