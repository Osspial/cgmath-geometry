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

pub mod iter;

use {AbsDistance, MulDiv, BaseScalarGeom, BaseVectorGeom, Dimensionality, D1, D2, D3};
use cgmath::*;

use line::{Linear, Segment};

use std::ops::{Add, Sub};
use num_traits::{Bounded, NumCast, ToPrimitive};

use self::iter::NearestPointsIter;

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct DimsBox<D: Dimensionality<S>, S: BaseScalarGeom> {
    pub dims: D::Vector
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct OffsetBox<D: Dimensionality<S>, S: BaseScalarGeom> {
    pub origin: D::Point,
    pub dims: D::Vector
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct BoundBox<D: Dimensionality<S>, S: BaseScalarGeom> {
    pub min: D::Point,
    pub max: D::Point
}


pub trait GeoBox: Sized
    where <Self::D as Dimensionality<Self::Scalar>>::Vector: VectorSpace<Scalar=Self::Scalar>,
          <Self::D as Dimensionality<Self::Scalar>>::Point: EuclideanSpace<Scalar=Self::Scalar, Diff=<Self::D as Dimensionality<Self::Scalar>>::Vector>
{
    type Scalar: BaseScalarGeom;
    type D: Dimensionality<Self::Scalar>;

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
    fn dims(&self) -> DimsBox<Self::D, Self::Scalar> {
        DimsBox::new(self.max() - self.min())
    }

    #[inline]
    fn width(&self) -> Self::Scalar {
        match d!(Point::len()) >= 1 {
            true => self.dims().dims[0],
            false => Self::Scalar::zero()
        }
    }
    #[inline]
    fn height(&self) -> Self::Scalar {
        match d!(Point::len()) >= 2 {
            true => self.dims().dims[1],
            false => Self::Scalar::zero()
        }
    }
    #[inline]
    fn depth(&self) -> Self::Scalar {
        match d!(Point::len()) >= 3 {
            true => self.dims().dims[2],
            false => Self::Scalar::zero()
        }
    }

    #[inline]
    fn center(&self) -> d!(Point) {
        self.min() + self.dims().dims / (Self::Scalar::one() + Self::Scalar::one())
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
        where L: Linear<Scalar=Self::Scalar, D=Self::D>,
              Self::Scalar: Bounded
    {
        struct TS;
        trait TypeSwitch<S>
            where S: BaseScalarGeom
        {
            fn intersect_ts<D, R, L>(rect: &R, line: L) -> (Option<D::Point>, Option<D::Point>)
                    where D: Dimensionality<S>,
                          D::Vector: BaseVectorGeom<D=D>,// VectorSpace<Scalar=S>,
                          D::Point: EuclideanSpace<Scalar=S, Diff=D::Vector>,
                          R: GeoBox<Scalar=S, D=D> + ?Sized,
                          L: Linear<Scalar=S, D=D>;
        }

        impl<S> TypeSwitch<S> for TS
            where S: BaseScalarGeom
        {
            default_if_nightly!{
                fn intersect_ts<D, R, L>(rect: &R, line: L) -> (Option<D::Point>, Option<D::Point>)
                        where D: Dimensionality<S>,
                            D::Vector: BaseVectorGeom<D=D>,
                            D::Point: EuclideanSpace<Scalar=S, Diff=D::Vector>,
                            R: GeoBox<Scalar=S, D=D> + ?Sized,
                            L: Linear<Scalar=S, D=D>
                {
                    let zero = S::zero();
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
        }

        // Optimized version for floating-point numbers
        #[cfg(feature="nightly")]
        impl<S> TypeSwitch<S> for TS
            where S: BaseFloat + BaseScalarGeom
        {
            fn intersect_ts<D, R, L>(rect: &R, line: L) -> (Option<D::Point>, Option<D::Point>)
                    where D: Dimensionality<S>,
                          D::Vector: BaseVectorGeom<D=D>,
                          D::Point: EuclideanSpace<Scalar=S, Diff=D::Vector>,
                          R: GeoBox<Scalar=S, D=D> + ?Sized,
                          L: Linear<Scalar=S, D=D>
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

        TS::intersect_ts(self, line)
    }

    fn clamp(&self, mut point: d!(Point)) -> d!(Point) {
        let min = self.min();
        let max = self.max();
        for i in 0..d!(Point::len()) {
            point[i] = num_traits::clamp(point[i], min[i], max[i]);
        }
        point
    }

    fn nearest_points(self, point: d!(Point)) -> NearestPointsIter<Self> {
        match self.contains(point) {
            false => NearestPointsIter {
                point: self.clamp(point),
                b: self,
                closest_sides: !0,
            },
            true => {
                let min = self.min();
                let max = self.max();

                let mut closest_axis_dist = <Self::Scalar as AbsDistance>::Abs::max_value();

                let mut closest_sides = 0;

                for i in 0..d!(Point::len()) {
                    let axis_dist_from_min = min[i].abs_distance(point[i]);
                    let axis_dist_from_max = max[i].abs_distance(point[i]);
                    let axis_dist = match axis_dist_from_min < axis_dist_from_max {
                        true => axis_dist_from_min,
                        false => axis_dist_from_max
                    };

                    if axis_dist < closest_axis_dist {
                        closest_sides = 0;
                        closest_axis_dist = axis_dist;
                    }

                    if axis_dist <= closest_axis_dist {
                        if axis_dist_from_min == axis_dist_from_max {
                            closest_sides |= 0b11 << (i * 2);
                        } else if axis_dist_from_min < axis_dist_from_max {
                            closest_sides |= 0b01 << (i * 2);
                        } else if axis_dist_from_max < axis_dist_from_min {
                            closest_sides |= 0b10 << (i * 2);
                        }
                    }
                }

                NearestPointsIter {
                    b: self,
                    point,
                    closest_sides,
                }
            }
        }
    }
}

macro_rules! inherent_impl_dims_offset {
    ($D:ident; $new:ident; $($origin:ident, $dim:ident),+) => {
        impl<S: BaseScalarGeom> DimsBox<$D, S> {
            #[inline]
            pub fn $new($($dim: S),+) -> DimsBox<$D, S> {
                DimsBox {
                    dims: <$D as Dimensionality<S>>::Vector::new($($dim),+)
                }
            }

            #[inline]
            pub fn cast<T>(&self) -> Option<DimsBox<$D, T>>
                where T: NumCast + BaseScalarGeom,
                      S: ToPrimitive
            {
                Some(DimsBox {
                    dims: <$D as Dimensionality<S>>::Vector::cast(&self.dims)?
                })
            }
        }

        impl<S: BaseScalarGeom> OffsetBox<$D, S> {
            #[inline]
            pub fn $new($($origin: S,)+ $($dim: S),+) -> OffsetBox<$D, S> {
                OffsetBox {
                    origin: <$D as Dimensionality<S>>::Point::new($($origin),+),
                    dims: <$D as Dimensionality<S>>::Vector::new($($dim),+)
                }
            }

            #[inline]
            pub fn cast<T>(&self) -> Option<OffsetBox<$D, T>>
                where T: NumCast + BaseScalarGeom,
                      S: ToPrimitive
            {
                Some(OffsetBox {
                    origin: <$D as Dimensionality<S>>::Point::cast(&self.origin)?,
                    dims: <$D as Dimensionality<S>>::Vector::cast(&self.dims)?
                })
            }
        }
    }
}

macro_rules! inherent_impl_bounds {
    ($D:ident; $new:ident; ($($min:ident),+), ($($max:ident),+)) => {
        impl<S: BaseScalarGeom> BoundBox<$D, S> {
            #[inline]
            pub fn $new($($min: S),+, $($max: S),+) -> BoundBox<$D, S> {
                BoundBox {
                    min: <$D as Dimensionality<S>>::Point::new($($min),+),
                    max: <$D as Dimensionality<S>>::Point::new($($max),+)
                }
            }

            #[inline]
            pub fn cast<T>(&self) -> Option<BoundBox<$D, T>>
                where T: NumCast + BaseScalarGeom,
                      S: ToPrimitive
            {
                Some(BoundBox {
                    min: <$D as Dimensionality<S>>::Point::cast(&self.min)?,
                    max: <$D as Dimensionality<S>>::Point::cast(&self.max)?,
                })
            }
        }
    }
}

impl<D: Dimensionality<S>, S: BaseScalarGeom> DimsBox<D, S> {
    #[inline]
    pub fn new(dims: D::Vector) -> DimsBox<D, S> {
        DimsBox{ dims }
    }
}
impl<D: Dimensionality<S>, S: BaseScalarGeom> OffsetBox<D, S> {
    #[inline]
    pub fn new(origin: D::Point, dims: D::Vector) -> OffsetBox<D, S> {
        OffsetBox{ origin, dims }
    }
}
impl<D: Dimensionality<S>, S: BaseScalarGeom> BoundBox<D, S> {
    #[inline]
    pub fn new(min: D::Point, max: D::Point) -> BoundBox<D, S> {
        BoundBox{ min, max }
    }
}

inherent_impl_dims_offset!(D1; new1; origin_x, width);
inherent_impl_dims_offset!(D2; new2; origin_x, width, origin_y, height);
inherent_impl_dims_offset!(D3; new3; origin_x, width, origin_y, height, origin_z, depth);

inherent_impl_bounds!(D1; new1; (min_x), (max_x));
inherent_impl_bounds!(D2; new2; (min_x, min_y), (max_x, max_y));
inherent_impl_bounds!(D3; new3; (min_x, min_y, min_z), (max_x, max_y, max_z));

impl<D, S> GeoBox for DimsBox<D, S>
    where S: BaseScalarGeom,
          D: Dimensionality<S>,
          D::Vector: VectorSpace<Scalar=S>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    type Scalar = S;
    type D = D;

    #[inline]
    fn from_bounds(min: D::Point, max: D::Point) -> DimsBox<D, S> {
        DimsBox {
            dims: max - min
        }
    }

    #[inline]
    fn min(&self) -> D::Point {D::Point::from_value(S::zero())}
    #[inline]
    fn dims(&self) -> DimsBox<D, S> {DimsBox{dims: self.dims}}
}

impl<D, S> Bounded for DimsBox<D, S>
    where S: BaseScalarGeom,
          D: Dimensionality<S>,
          D::Vector: Bounded
{
    #[inline]
    fn min_value() -> DimsBox<D, S> {
        DimsBox {
            dims: D::Vector::min_value()
        }
    }

    #[inline]
    fn max_value() -> DimsBox<D, S> {
        DimsBox {
            dims: D::Vector::max_value()
        }
    }
}

impl<D, S> GeoBox for OffsetBox<D, S>
    where S: BaseScalarGeom,
          D: Dimensionality<S>,
          D::Vector: VectorSpace<Scalar=S>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    type Scalar = S;
    type D = D;

    #[inline]
    fn from_bounds(min: D::Point, max: D::Point) -> OffsetBox<D, S> {
        OffsetBox {
            origin: min,
            dims: max - min
        }
    }

    #[inline]
    fn min(&self) -> D::Point {self.origin}
    #[inline]
    fn dims(&self) -> DimsBox<D, S> {DimsBox::new(self.dims)}
}

impl<D, S> GeoBox for BoundBox<D, S>
    where S: BaseScalarGeom,
          D: Dimensionality<S>,
          D::Vector: VectorSpace<Scalar=S>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    type Scalar = S;
    type D = D;

    #[inline]
    fn from_bounds(min: D::Point, max: D::Point) -> BoundBox<D, S> {
        BoundBox {
            min, max
        }
    }

    #[inline]
    fn min(&self) -> D::Point {self.min}
    #[inline]
    fn max(&self) -> D::Point {self.max}
}

impl<D, S> From<DimsBox<D, S>> for OffsetBox<D, S>
    where S: BaseScalarGeom,
          D: Dimensionality<S>
{
    #[inline]
    fn from(rect: DimsBox<D, S>) -> OffsetBox<D, S> {
        OffsetBox {
            origin: D::Point::from_value(S::zero()),
            dims: rect.dims
        }
    }
}

impl<D, S> From<DimsBox<D, S>> for BoundBox<D, S>
    where S: BaseScalarGeom,
          D: Dimensionality<S>,
          D::Vector: VectorSpace<Scalar=S>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    #[inline]
    fn from(rect: DimsBox<D, S>) -> BoundBox<D, S> {
        BoundBox {
            min: D::Point::from_value(S::zero()),
            max: D::Point::from_vec(rect.dims)
        }
    }
}

impl<D, S> From<OffsetBox<D, S>> for BoundBox<D, S>
    where S: BaseScalarGeom,
          D: Dimensionality<S>,
          D::Vector: VectorSpace<Scalar=S>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    #[inline]
    fn from(rect: OffsetBox<D, S>) -> BoundBox<D, S> {
        BoundBox {
            min: rect.origin,
            max: D::Point::add(rect.origin, rect.dims)
        }
    }
}

impl<D, S> From<BoundBox<D, S>> for OffsetBox<D, S>
    where S: BaseScalarGeom,
          D: Dimensionality<S>,
          D::Vector: VectorSpace<Scalar=S>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    #[inline]
    fn from(rect: BoundBox<D, S>) -> OffsetBox<D, S> {
        OffsetBox {
            origin: rect.min,
            dims: D::Point::to_vec(rect.max) - D::Point::to_vec(rect.min)
        }
    }
}

impl<D, S> Add<D::Vector> for OffsetBox<D, S>
    where S: BaseScalarGeom,
          D: Dimensionality<S>,
          D::Vector: VectorSpace<Scalar=S>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    type Output = Self;
    #[inline]
    fn add(mut self, rhs: D::Vector) -> OffsetBox<D, S> {
        self.origin = self.origin + rhs;
        self
    }
}

impl<D, S> Sub<D::Vector> for OffsetBox<D, S>
    where S: BaseScalarGeom,
          D: Dimensionality<S>,
          D::Vector: VectorSpace<Scalar=S>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    type Output = Self;
    #[inline]
    fn sub(mut self, rhs: D::Vector) -> OffsetBox<D, S> {
        self.origin = self.origin - rhs;
        self
    }
}

impl<D, S> Add<D::Vector> for BoundBox<D, S>
    where S: BaseScalarGeom,
          D: Dimensionality<S>,
          D::Vector: VectorSpace<Scalar=S>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    type Output = Self;
    #[inline]
    fn add(mut self, rhs: D::Vector) -> BoundBox<D, S> {
        self.min = self.min + rhs;
        self.max = self.max + rhs;
        self
    }
}

impl<D, S> Sub<D::Vector> for BoundBox<D, S>
    where S: BaseScalarGeom,
          D: Dimensionality<S>,
          D::Vector: VectorSpace<Scalar=S>,
          D::Point: EuclideanSpace<Diff=D::Vector>
{
    type Output = Self;
    #[inline]
    fn sub(mut self, rhs: D::Vector) -> BoundBox<D, S> {
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

    #[test]
    fn contains() {
        let rect = BoundBox::new2(30, 10, 50, 20);

        // Test borders
        assert!(rect.contains(Point2::new(30, 10)));
        assert!(rect.contains(Point2::new(40, 10)));
        assert!(rect.contains(Point2::new(50, 10)));
        assert!(rect.contains(Point2::new(50, 15)));
        assert!(rect.contains(Point2::new(50, 20)));
        assert!(rect.contains(Point2::new(40, 20)));
        assert!(rect.contains(Point2::new(30, 20)));
        assert!(rect.contains(Point2::new(30, 15)));

        assert!(rect.contains(Point2::new(40, 15)));

        assert!(!rect.contains(Point2::new(25, 5)));
        assert!(!rect.contains(Point2::new(40, 5)));
        assert!(!rect.contains(Point2::new(55, 5)));
        assert!(!rect.contains(Point2::new(55, 15)));
        assert!(!rect.contains(Point2::new(55, 25)));
        assert!(!rect.contains(Point2::new(40, 25)));
        assert!(!rect.contains(Point2::new(25, 25)));
        assert!(!rect.contains(Point2::new(25, 15)));
    }

    #[test]
    fn test_nearest_points() {
        let rect = BoundBox::new2(30, 10, 50, 20);

        // Outside points
        assert_eq!(
            rect.nearest_points(Point2::new(0, 0)).collect::<Vec<_>>(),
            &[Point2::new(30, 10)]
        );
        assert_eq!(
            rect.nearest_points(Point2::new(0, 20)).collect::<Vec<_>>(),
            &[Point2::new(30, 20)]
        );

        // Inside points
        assert_eq!(
            rect.nearest_points(Point2::new(33, 15)).collect::<Vec<_>>(),
            &[Point2::new(30, 15)]
        );
        assert_eq!(
            rect.nearest_points(Point2::new(47, 15)).collect::<Vec<_>>(),
            &[Point2::new(50, 15)]
        );
        assert_eq!(
            rect.nearest_points(Point2::new(40, 16)).collect::<Vec<_>>(),
            &[Point2::new(40, 20)]
        );
        assert_eq!(
            rect.nearest_points(Point2::new(40, 14)).collect::<Vec<_>>(),
            &[Point2::new(40, 10)]
        );
        assert_eq!(
            rect.nearest_points(Point2::new(35, 15)).collect::<Vec<_>>(),
            &[Point2::new(30, 15), Point2::new(35, 10), Point2::new(35, 20)]
        );
        assert_eq!(
            rect.nearest_points(Point2::new(40, 15)).collect::<Vec<_>>(),
            &[Point2::new(40, 10), Point2::new(40, 20)]
        );
        assert_eq!(
            rect.nearest_points(Point2::new(45, 15)).collect::<Vec<_>>(),
            &[Point2::new(50, 15), Point2::new(45, 10), Point2::new(45, 20)]
        );
    }
}
