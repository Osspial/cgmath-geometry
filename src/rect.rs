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

use {MulDiv, BaseScalarGeom, AbsDistance, BasePointGeom, BaseVectorGeom, Dimensionality, D1, D2, D3};
use cgmath::*;

use line::{Linear, Segment};

use std::ops::{Add, Sub};
use num_traits::{Bounded, NumCast, ToPrimitive};

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct DimsBox<D: Dimensionality> {
    pub dims: D::Vector
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct OffsetBox<D: Dimensionality> {
    pub origin: D::Point,
    pub dims: D::Vector
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct BoundBox<D: Dimensionality> {
    pub min: D::Point,
    pub max: D::Point
}

pub trait GeoBox
    where <Self::D as Dimensionality>::Vector: VectorSpace<Scalar=<Self::D as Dimensionality>::Scalar>,
          <Self::D as Dimensionality>::Point: EuclideanSpace<Diff=<Self::D as Dimensionality>::Vector>
{
    type D: Dimensionality;

    fn from_bounds(min: d!(Point), max: d!(Point)) -> Self;

    #[inline]
    fn min(&self) -> d!(Point) {
        self.max() - self.dims().dims
    }

    #[inline]
    fn max(&self) -> d!(Point) {
        d!(Point::add(self.min(), self.dims().dims))
    }

    #[inline]
    fn dims(&self) -> DimsBox<Self::D> {
        DimsBox::new(self.max() - self.min())
    }

    #[inline]
    fn width(&self) -> d!(Scalar) {
        match d!(Point::len()) >= 1 {
            true => self.dims().dims[0],
            false => d!(Scalar::zero())
        }
    }
    #[inline]
    fn height(&self) -> d!(Scalar) {
        match d!(Point::len()) >= 2 {
            true => self.dims().dims[1],
            false => d!(Scalar::zero())
        }
    }
    #[inline]
    fn depth(&self) -> d!(Scalar) {
        match d!(Point::len()) >= 3 {
            true => self.dims().dims[2],
            false => d!(Scalar::zero())
        }
    }

    #[inline]
    fn center(&self) -> d!(Point) {
        self.min() + self.dims().dims / (d!(Scalar::one()) + d!(Scalar::one()))
    }

    #[inline]
    fn contains(&self, point: d!(Point)) -> bool {
        let min = self.min();
        let max = self.max();

        let mut contains = true;
        for i in 0..d!(Point::len()) {
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

        for i in 0..d!(Point::len()) {
            min[i] = ::cmp_max(s_min[i], o_min[i]);
            max[i] = ::cmp_min(s_max[i], o_max[i]);

            if max[i] < min[i] {
                return None;
            }
        }

        Some(Self::from_bounds(min, max))
    }

    fn intersect_line<L>(&self, line: L) -> (Option<d!(Point)>, Option<d!(Point)>)
        where L: Linear<D=Self::D>,
              d!(Scalar): Bounded
    {
        struct TS;
        trait TypeSwitch<S>
            where S: BaseScalarGeom
        {
            fn intersect_ts<D, R, L>(rect: &R, line: L) -> (Option<D::Point>, Option<D::Point>)
                    where D: Dimensionality<Scalar=S>,
                          D::Vector: BaseVectorGeom<D=D>,// VectorSpace<Scalar=D::Scalar>,
                          D::Point: BasePointGeom<D=D>,// EuclideanSpace<Diff=D::Vector>,
                          R: GeoBox<D=D> + ?Sized,
                          L: Linear<D=D>;
        }

        impl<S> TypeSwitch<S> for TS
            where S: BaseFloat + BaseScalarGeom
        {
            fn intersect_ts<D, R, L>(rect: &R, line: L) -> (Option<D::Point>, Option<D::Point>)
                    where D: Dimensionality<Scalar=S>,
                          D::Vector: VectorSpace<Scalar=D::Scalar>,
                          D::Point: EuclideanSpace<Diff=D::Vector>,
                          R: GeoBox<D=D> + ?Sized,
                          L: Linear<D=D>
            {
                let line_origin = line.origin();
                let dir = line.dir();
                let dir_recip = line.dir_recip();
                let (rect_min, rect_max) = (rect.min(), rect.max());

                let (mut t_min, mut t_max) = (S::neg_infinity(), S::infinity());

                for i in 0..D::Point::len() {
                    let t_enter = (rect_min[i] - line_origin[i]) * dir_recip[i];
                    let t_exit = (rect_max[i] - line_origin[i]) * dir_recip[i];
                    t_min = t_min.max(t_enter.min(t_exit));
                    t_max = t_max.min(t_enter.max(t_exit));
                }

                if t_max < t_min {
                    (None, None)
                } else {
                    let t_of_point = |point: D::Point| {
                        let mut t = S::zero();
                        for i in 0..D::Point::len() {
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
            }
        }

        impl<S> TypeSwitch<S> for TS
            where S: BaseScalarGeom
        {
            default fn intersect_ts<D, R, L>(rect: &R, line: L) -> (Option<D::Point>, Option<D::Point>)
                    where D: Dimensionality<Scalar=S>,
                          D::Vector: VectorSpace<Scalar=D::Scalar>,
                          D::Point: EuclideanSpace<Diff=D::Vector>,
                          R: GeoBox<D=D> + ?Sized,
                          L: Linear<D=D>
            {
                let zero = D::Scalar::zero();
                let min = rect.min();
                let max = rect.max();

                let Segment{ start, end } = line.clip_to_scalar_bounds();
                let (mut enter, mut exit) = (start, end);
                let (mut enter_valid, mut exit_valid) = (false, false);
                let dir = line.dir();

                for i in 0..D::Point::len() {
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

                for i in 0..D::Point::len() {
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
            }
        }

        TS::intersect_ts(self, line)
    }
}

macro_rules! inherent_impl_dims_offset {
    ($D:ident; $new:ident; $($origin:ident, $dim:ident),+) => {
        impl<S: BaseScalarGeom> DimsBox<$D<S>> {
            #[inline]
            pub fn $new($($dim: S),+) -> DimsBox<$D<S>> {
                DimsBox {
                    dims: <$D<S> as Dimensionality>::Vector::new($($dim),+)
                }
            }

            #[inline]
            pub fn cast<T>(&self) -> Option<DimsBox<$D<T>>>
                where T: NumCast + BaseScalarGeom,
                      S: ToPrimitive
            {
                Some(DimsBox {
                    dims: <$D<S> as Dimensionality>::Vector::cast(&self.dims)?
                })
            }
        }

        impl<S: BaseScalarGeom> OffsetBox<$D<S>> {
            #[inline]
            pub fn $new($($origin: S,)+ $($dim: S),+) -> OffsetBox<$D<S>> {
                OffsetBox {
                    origin: <$D<S> as Dimensionality>::Point::new($($origin),+),
                    dims: <$D<S> as Dimensionality>::Vector::new($($dim),+)
                }
            }

            #[inline]
            pub fn cast<T>(&self) -> Option<OffsetBox<$D<T>>>
                where T: NumCast + BaseScalarGeom,
                      S: ToPrimitive
            {
                Some(OffsetBox {
                    origin: <$D<S> as Dimensionality>::Point::cast(&self.origin)?,
                    dims: <$D<S> as Dimensionality>::Vector::cast(&self.dims)?
                })
            }
        }
    }
}

macro_rules! inherent_impl_bounds {
    ($D:ident; $new:ident; ($($min:ident),+), ($($max:ident),+)) => {
        impl<S: BaseScalarGeom> BoundBox<$D<S>> {
            #[inline]
            pub fn $new($($min: S),+, $($max: S),+) -> BoundBox<$D<S>> {
                BoundBox {
                    min: <$D<S> as Dimensionality>::Point::new($($min),+),
                    max: <$D<S> as Dimensionality>::Point::new($($max),+)
                }
            }

            #[inline]
            pub fn cast<T>(&self) -> Option<BoundBox<$D<T>>>
                where T: NumCast + BaseScalarGeom,
                      S: ToPrimitive
            {
                Some(BoundBox {
                    min: <$D<S> as Dimensionality>::Point::cast(&self.min)?,
                    max: <$D<S> as Dimensionality>::Point::cast(&self.max)?,
                })
            }
        }
    }
}

impl<D: Dimensionality> DimsBox<D> {
    #[inline]
    pub fn new(dims: D::Vector) -> DimsBox<D> {
        DimsBox{ dims }
    }
}
impl<D: Dimensionality> OffsetBox<D> {
    #[inline]
    pub fn new(origin: D::Point, dims: D::Vector) -> OffsetBox<D> {
        OffsetBox{ origin, dims }
    }
}
impl<D: Dimensionality> BoundBox<D> {
    #[inline]
    pub fn new(min: D::Point, max: D::Point) -> BoundBox<D> {
        BoundBox{ min, max }
    }
}

inherent_impl_dims_offset!(D1; new1; origin_x, width);
inherent_impl_dims_offset!(D2; new2; origin_x, width, origin_y, height);
inherent_impl_dims_offset!(D3; new3; origin_x, width, origin_y, height, origin_z, depth);

inherent_impl_bounds!(D1; new1; (min_x), (max_x));
inherent_impl_bounds!(D2; new2; (min_x, min_y), (max_x, max_y));
inherent_impl_bounds!(D3; new3; (min_x, min_y, min_z), (max_x, max_y, max_z));

impl<D> GeoBox for DimsBox<D>
    where D: Dimensionality,
          D::Vector: VectorSpace<Scalar=D::Scalar>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    type D = D;

    #[inline]
    fn from_bounds(min: D::Point, max: D::Point) -> DimsBox<D> {
        DimsBox {
            dims: max - min
        }
    }

    #[inline]
    fn min(&self) -> D::Point {D::Point::from_value(D::Scalar::zero())}
    #[inline]
    fn dims(&self) -> DimsBox<D> {DimsBox{dims: self.dims}}
}

impl<D> Bounded for DimsBox<D>
    where D: Dimensionality,
          D::Vector: Bounded
{
    #[inline]
    fn min_value() -> DimsBox<D> {
        DimsBox {
            dims: D::Vector::min_value()
        }
    }

    #[inline]
    fn max_value() -> DimsBox<D> {
        DimsBox {
            dims: D::Vector::max_value()
        }
    }
}

impl<D> GeoBox for OffsetBox<D>
    where D: Dimensionality,
          D::Vector: VectorSpace<Scalar=D::Scalar>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    type D = D;

    #[inline]
    fn from_bounds(min: D::Point, max: D::Point) -> OffsetBox<D> {
        OffsetBox {
            origin: min,
            dims: max - min
        }
    }

    #[inline]
    fn min(&self) -> D::Point {self.origin}
    #[inline]
    fn dims(&self) -> DimsBox<D> {DimsBox::new(self.dims)}
}

impl<D> GeoBox for BoundBox<D>
    where D: Dimensionality,
          D::Vector: VectorSpace<Scalar=D::Scalar>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    type D = D;

    #[inline]
    fn from_bounds(min: D::Point, max: D::Point) -> BoundBox<D> {
        BoundBox {
            min, max
        }
    }

    #[inline]
    fn min(&self) -> D::Point {self.min}
    #[inline]
    fn max(&self) -> D::Point {self.max}
}

impl<D> From<DimsBox<D>> for OffsetBox<D>
    where D: Dimensionality
{
    #[inline]
    fn from(rect: DimsBox<D>) -> OffsetBox<D> {
        OffsetBox {
            origin: D::Point::from_value(D::Scalar::zero()),
            dims: rect.dims
        }
    }
}

impl<D> From<DimsBox<D>> for BoundBox<D>
    where D: Dimensionality,
          D::Vector: VectorSpace<Scalar=D::Scalar>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    #[inline]
    fn from(rect: DimsBox<D>) -> BoundBox<D> {
        BoundBox {
            min: D::Point::from_value(D::Scalar::zero()),
            max: D::Point::from_vec(rect.dims)
        }
    }
}

impl<D> From<OffsetBox<D>> for BoundBox<D>
    where D: Dimensionality,
          D::Vector: VectorSpace<Scalar=D::Scalar>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    #[inline]
    fn from(rect: OffsetBox<D>) -> BoundBox<D> {
        BoundBox {
            min: rect.origin,
            max: D::Point::add(rect.origin, rect.dims)
        }
    }
}

impl<D> From<BoundBox<D>> for OffsetBox<D>
    where D: Dimensionality,
          D::Vector: VectorSpace<Scalar=D::Scalar>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    #[inline]
    fn from(rect: BoundBox<D>) -> OffsetBox<D> {
        OffsetBox {
            origin: rect.min,
            dims: D::Point::to_vec(rect.max) - D::Point::to_vec(rect.min)
        }
    }
}

impl<D> Add<D::Vector> for OffsetBox<D>
    where D: Dimensionality,
          D::Vector: VectorSpace<Scalar=D::Scalar>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    type Output = Self;
    #[inline]
    fn add(mut self, rhs: D::Vector) -> OffsetBox<D> {
        self.origin = self.origin + rhs;
        self
    }
}

impl<D> Sub<D::Vector> for OffsetBox<D>
    where D: Dimensionality,
          D::Vector: VectorSpace<Scalar=D::Scalar>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    type Output = Self;
    #[inline]
    fn sub(mut self, rhs: D::Vector) -> OffsetBox<D> {
        self.origin = self.origin - rhs;
        self
    }
}

impl<D> Add<D::Vector> for BoundBox<D>
    where D: Dimensionality,
          D::Vector: VectorSpace<Scalar=D::Scalar>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    type Output = Self;
    #[inline]
    fn add(mut self, rhs: D::Vector) -> BoundBox<D> {
        self.min = self.min + rhs;
        self.max = self.max + rhs;
        self
    }
}

impl<D> Sub<D::Vector> for BoundBox<D>
    where D: Dimensionality,
          D::Vector: VectorSpace<Scalar=D::Scalar>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    type Output = Self;
    #[inline]
    fn sub(mut self, rhs: D::Vector) -> BoundBox<D> {
        self.min = D::Point::sub(self.min, rhs);
        self.max = D::Point::sub(self.max, rhs);
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
