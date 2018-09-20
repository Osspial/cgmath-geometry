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

use {MulDiv, BaseScalarGeom, Intersect, Intersection, AbsDistance, Dimensionality, D1, D2, D3};
use cgmath::*;
use num_traits::{Bounded, Float, Signed};
use std::cmp::{Ordering, PartialEq, Eq};
use std::ops::Mul;
use rect::{BoundBox, GeoBox};

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct Ray<D: Dimensionality> {
    pub origin: D::Point,
    pub dir: D::Vector
}

// Eq and PartialEq manually implemented
#[repr(C)]
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct Segment<D: Dimensionality> {
    pub start: D::Point,
    pub end: D::Point
}

// Eq and PartialEq manually implemented
#[repr(C)]
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct Line<D: Dimensionality> {
    pub origin: D::Point,
    pub dir: D::Vector
}

pub trait Linear {
    type D: Dimensionality;

    fn origin(&self) -> d!(Point);
    fn dir(&self) -> d!(Vector);
    #[inline]
    fn dir_recip(&self) -> d!(Vector)
        where d!(Scalar): Float
    {
        let mut d = self.dir();
        for i in 0..d!(Vector::len()) {
            d[i] = d[i].recip();
        }
        d
    }
    fn start(&self) -> Option<d!(Point)>;
    fn end(&self) -> Option<d!(Point)>;

    fn bounding_box(&self) -> BoundBox<Self::D> {
        let l = self.start().expect("Manual bounding box impl required");
        let r = self.end().expect("Manual bounding box impl required");

        let mut min = d!(Point::from_value(d!(Scalar::zero())));
        let mut max = min;
        for i in 0..d!(Point::len()) {
            min[i] = ::cmp_min(l[i], r[i]);
            max[i] = ::cmp_max(l[i], r[i]);
        }

        BoundBox::from_bounds(min, max)
    }

    fn clip_to_scalar_bounds(&self) -> Segment<Self::D> {
        let origin = self.origin();
        let dir = self.dir();
        let zero = d!(Scalar::zero());

        // to_end is whether to project to the end of the line (along the direction vector) or to
        // the start of the line (against the direction vector). Changing that changes some
        // comparisons.
        let project = |to_end: bool| {
            let mut axis_closest = usize::max_value();
            let mut axis_closest_threshold = d!(Scalar::max_value().to_abs());
            let mut axis_closest_muldiv = (zero.to_abs(), zero.to_abs());

            for i in 0..d!(Point::len()) {
                let dir_cmp = match to_end {
                    true => dir[i].partial_cmp(&zero),
                    false => zero.partial_cmp(&dir[i])
                };
                let axis_dist = match dir_cmp {
                    Some(Ordering::Less) => origin[i].abs_distance(d!(Scalar::min_value())),
                    Some(Ordering::Greater) => origin[i].abs_distance(d!(Scalar::max_value())),
                    Some(Ordering::Equal) |
                    None                 => d!(Scalar::max_value().to_abs())
                };

                let axis_threshold = axis_dist / dir[i].to_abs();
                if axis_threshold <= axis_closest_threshold {
                    axis_closest = i;
                    axis_closest_threshold = axis_threshold;
                    axis_closest_muldiv = (axis_dist, dir[i].to_abs());
                }
            }

            if axis_closest == usize::max_value() {
                panic!("direction equal to zero");
            } else {
                let mut origin_projected = origin;
                for i in 0..d!(Point::len()) {
                    let dir_mul_div = dir[i].to_abs().mul_div(axis_closest_muldiv.0, axis_closest_muldiv.1);
                    // `... ^ !to_end` flips the sign if we're going to the start (switch add to sub, and sub to add)
                    origin_projected[i] = match (zero <= dir[i]) ^ !to_end {
                        true => origin[i].add_abs(dir_mul_div),
                        false => origin[i].sub_abs(dir_mul_div),
                    };
                }
                origin_projected
            }
        };

        Segment {
            start: self.start().unwrap_or_else(|| project(false)),
            end: self.end().unwrap_or_else(|| project(true))
        }
    }
}

impl<D> Linear for Ray<D>
    where D: Dimensionality
{
    type D = D;

    #[inline]
    fn origin(&self) -> D::Point {self.origin}
    #[inline]
    fn dir(&self) -> D::Vector {self.dir}
    #[inline]
    fn start(&self) -> Option<D::Point> {
        Some(self.origin)
    }
    #[inline]
    fn end(&self) -> Option<D::Point> {
        None
    }

    #[inline]
    default fn bounding_box(&self) -> BoundBox<D> {
        let mut outer_bound = d!(Point::from_vec(self.dir()));
        for i in 0..D::Point::len() {
            outer_bound[i] = match outer_bound[i].partial_cmp(&D::Scalar::zero()) {
                Some(Ordering::Less) => D::Scalar::min_value(),
                Some(Ordering::Greater) => D::Scalar::max_value(),
                None |
                Some(Ordering::Equal) => self.origin[i],
            };
        }

        let (l, r) = (self.origin, outer_bound);
        let mut min = D::Point::from_value(D::Scalar::zero());
        let mut max = min;
        for i in 0..D::Point::len() {
            min[i] = ::cmp_min(l[i], r[i]);
            max[i] = ::cmp_max(l[i], r[i]);
        }

        BoundBox::from_bounds(min, max)
    }
}

impl<D> Linear for Ray<D>
    where D: Dimensionality,
          D::Scalar: BaseFloat
{
    #[inline]
    fn bounding_box(&self) -> BoundBox<D> {
        let outer_bound = d!(Point::from_vec(self.dir()).mul_element_wise(D::Scalar::infinity()));

        let (l, r) = (self.origin, outer_bound);
        let mut min = D::Point::from_value(D::Scalar::zero());
        let mut max = min;
        for i in 0..D::Point::len() {
            min[i] = ::cmp_min(l[i], r[i]);
            max[i] = ::cmp_max(l[i], r[i]);
        }

        BoundBox::from_bounds(min, max)
    }
}

impl<D> Linear for Segment<D>
    where D: Dimensionality
{
    type D = D;

    #[inline]
    fn origin(&self) -> D::Point {
        self.start
    }
    #[inline]
    fn dir(&self) -> D::Vector {
        self.end - self.start
    }
    #[inline]
    fn start(&self) -> Option<D::Point> {
        Some(self.start)
    }
    #[inline]
    fn end(&self) -> Option<D::Point> {
        Some(self.end)
    }
}

impl<D> Linear for Line<D>
    where D: Dimensionality
{
    type D = D;

    fn origin(&self) -> D::Point {
        self.origin
    }

    fn dir(&self) -> D::Vector {
        self.dir
    }

    fn start(&self) -> Option<D::Point> {None}
    fn end(&self) -> Option<D::Point> {None}
}

impl<D: Dimensionality> Segment<D> {
    #[inline]
    pub fn new(start: D::Point, end: D::Point) -> Segment<D> {
        Segment{ start, end }
    }
}

macro_rules! inherent_impl_segment {
    ($D:ident; $new:ident; ($($start:ident),+), ($($end:ident),+)) => {
        impl<S: BaseScalarGeom> Segment<$D<S>> {
            #[inline]
            pub fn $new($($start: S),+, $($end: S),+) -> Segment<$D<S>> {
                Segment {
                    start: <$D<S> as Dimensionality>::Point::new($($start),+),
                    end: <$D<S> as Dimensionality>::Point::new($($end),+)
                }
            }
        }
    }
}
macro_rules! inherent_impl_ray {
    ($Name:ident; $D:ident; $new:ident; ($($origin:ident),+), ($($dir:ident),+)) => {
        impl<S: BaseScalarGeom> $Name<$D<S>> {
            #[inline]
            pub fn $new($($origin: S),+, $($dir: S),+) -> $Name<$D<S>> {
                $Name {
                    origin: <$D<S> as Dimensionality>::Point::new($($origin),+),
                    dir: <$D<S> as Dimensionality>::Vector::new($($dir),+)
                }
            }
        }
    }
}

inherent_impl_segment!(D1; new1; (start_x), (end_x));
inherent_impl_segment!(D2; new2; (start_x, start_y), (end_x, end_y));
inherent_impl_segment!(D3; new3; (start_x, start_y, start_z), (end_x, end_y, end_z));
inherent_impl_ray!(Ray; D1; new1; (origin_x), (dir_x));
inherent_impl_ray!(Ray; D2; new2; (origin_x, origin_y), (dir_x, dir_y));
inherent_impl_ray!(Ray; D3; new3; (origin_x, origin_y, origin_z), (dir_x, dir_y, dir_z));
inherent_impl_ray!(Line; D1; new1; (origin_x), (dir_x));
inherent_impl_ray!(Line; D2; new2; (origin_x, origin_y), (dir_x, dir_y));
inherent_impl_ray!(Line; D3; new3; (origin_x, origin_y, origin_z), (dir_x, dir_y, dir_z));

macro_rules! ld {
    ($($t:tt)*) => {<L::D as Dimensionality>::$($t)* };
}

impl<L, R> Intersect<R> for L
    where L: Linear,
          R: Linear<D=L::D>,
          L::D: Dimensionality<Point=Point2<ld!(Scalar)>, Vector=Vector2<ld!(Scalar)>>
{
    type Intersection = ld!(Point);
    fn intersect(self, rhs: R) -> Intersection<ld!(Point)> {
        let (lo, ro) = (self.origin(), rhs.origin());
        let (ld, rd) = (self.dir(), rhs.dir());
        let slope_diff = rd.y*ld.x - rd.x*ld.y;

        if slope_diff == ld!(Scalar::zero()) {
            return if (ro.y-lo.y) * (ro.x-lo.x) == ld.x * ld.y || ro == lo {
                Intersection::Eq
            } else {
                Intersection::None
            };
        }
        let t = (rd.x*(lo.y-ro.y) - rd.y*(lo.x-ro.x))/slope_diff;
        let intersection = lo + ld * t;

        match self.bounding_box().intersect_rect(rhs.bounding_box()).map(|r| r.contains(intersection)) {
            Some(true) => Intersection::Some(intersection),
            _ => Intersection::None
        }
    }
}

// impl<P: EuclideanSpace> Intersect for Line<P>
//     where P::Scalar: BaseScalarGeom
// {
//     type Intersection = Point2<P::Scalar>;
//     fn intersect(self, rhs: Line<P>) -> Intersection<Point2<P::Scalar>> {
//         let (lo, ro) = (self.origin(), rhs.origin());
//         let (ld, rd) = (self.dir(), rhs.dir());
//         let slope_diff = rd.y*ld.x - rd.x*ld.y;

//         if slope_diff == P::Scalar::zero() {
//             return if (ro.y-lo.y) * (ro.x-lo.x) == ld.x * ld.y || ro == lo {
//                 Intersection::Eq
//             } else {
//                 Intersection::None
//             };
//         }
//         let t = (rd.x*(lo.y-ro.y) - rd.y*(lo.x-ro.x))/slope_diff;

//         Intersection::Some(lo + ld * t)
//     }
// }

impl<D: Dimensionality> PartialEq for Line<D>
    where D::Vector: PartialEq + Mul + Array<Element=D::Scalar>,
          D::Scalar: Mul
{
    #[inline]
    default fn eq(&self, other: &Line<D>) -> bool {
        self.dir == other.dir && (other.origin - self.origin).product() == self.dir.product()
    }
}

impl<D: Dimensionality> PartialEq for Line<D>
    where D::Vector: PartialEq + Mul + Array<Element=D::Scalar>,
          D::Scalar: Mul + Signed
{
    #[inline]
    fn eq(&self, other: &Line<D>) -> bool {
        self.dir == other.dir && (other.origin - self.origin).product().abs() == self.dir.product().abs()
    }
}
impl<D> Eq for Line<D>
    where D: Dimensionality,
          D::Vector: PartialEq + Mul + Array<Element=D::Scalar>,
          D::Scalar: Mul + Eq {}

