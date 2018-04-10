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

use {MulDiv, BaseScalarGeom, Intersect, Intersection, AbsDistance};
use cgmath::*;
use num_traits::{Bounded, Float, Signed};
use std::cmp::{Ordering, PartialEq, Eq};
use std::ops::Mul;
use rect::{BoundBox, GeoBox};

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct Ray<P: EuclideanSpace> {
    pub origin: P,
    pub dir: P::Diff
}

// Eq and PartialEq manually implemented
#[repr(C)]
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct Segment<P: EuclideanSpace> {
    pub start: P,
    pub end: P
}

// Eq and PartialEq manually implemented
#[repr(C)]
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct Line<P: EuclideanSpace> {
    pub origin: P,
    pub dir: P::Diff
}

pub trait Linear {
    type Scalar: BaseScalarGeom;
    type Point: EuclideanSpace<Scalar=Self::Scalar, Diff=Self::Vector> + ElementWise<Self::Scalar> + MulDiv<Self::Scalar>;
    type Vector: VectorSpace<Scalar=Self::Scalar> + Array<Element=Self::Scalar> + MulDiv + MulDiv<Self::Scalar> + ElementWise;

    fn origin(&self) -> Self::Point;
    fn dir(&self) -> Self::Vector;
    #[inline]
    fn dir_recip(&self) -> Self::Vector
        where Self::Scalar: Float
    {
        let mut d = self.dir();
        for i in 0..Self::Vector::len() {
            d[i] = d[i].recip();
        }
        d
    }
    fn start(&self) -> Option<Self::Point>;
    fn end(&self) -> Option<Self::Point>;

    fn bounding_box(&self) -> BoundBox<Self::Point> {
        let l = self.start().expect("Manual bounding box impl required");
        let r = self.end().expect("Manual bounding box impl required");

        let mut min = Self::Point::from_value(Self::Scalar::zero());
        let mut max = min;
        for i in 0..Self::Point::len() {
            min[i] = ::cmp_min(l[i], r[i]);
            max[i] = ::cmp_max(l[i], r[i]);
        }

        BoundBox::from_bounds(min, max)
    }

    fn clip_to_scalar_bounds(&self) -> Segment<Self::Point> {
        let origin = self.origin();
        let dir = self.dir();
        let zero = Self::Scalar::zero();

        // to_end is whether to project to the end of the line (along the direction vector) or to
        // the start of the line (against the direction vector). Changing that changes some
        // comparisons.
        let project = |to_end: bool| {
            let mut axis_closest = usize::max_value();
            let mut axis_closest_threshold = Self::Scalar::max_value().to_abs();
            let mut axis_closest_muldiv = (zero.to_abs(), zero.to_abs());

            for i in 0..Self::Point::len() {
                let dir_cmp = match to_end {
                    true => dir[i].partial_cmp(&zero),
                    false => zero.partial_cmp(&dir[i])
                };
                let axis_dist = match dir_cmp {
                    Some(Ordering::Less) => origin[i].abs_distance(Self::Scalar::min_value()),
                    Some(Ordering::Greater) => origin[i].abs_distance(Self::Scalar::max_value()),
                    Some(Ordering::Equal) |
                    None                 => Self::Scalar::max_value().to_abs()
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
                for i in 0..Self::Point::len() {
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

impl<P> Linear for Ray<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise
{
    type Scalar = P::Scalar;
    type Point = P;
    type Vector = P::Diff;

    #[inline]
    fn origin(&self) -> P {self.origin}
    #[inline]
    fn dir(&self) -> P::Diff {self.dir}
    #[inline]
    fn start(&self) -> Option<P> {
        Some(self.origin)
    }
    #[inline]
    fn end(&self) -> Option<P> {
        None
    }

    #[inline]
    default fn bounding_box(&self) -> BoundBox<P> {
        let mut outer_bound = Self::Point::from_vec(self.dir());
        for i in 0..Self::Point::len() {
            outer_bound[i] = match outer_bound[i].partial_cmp(&Self::Scalar::zero()) {
                Some(Ordering::Less) => Self::Scalar::min_value(),
                Some(Ordering::Greater) => Self::Scalar::max_value(),
                None |
                Some(Ordering::Equal) => self.origin[i],
            };
        }

        let (l, r) = (self.origin, outer_bound);
        let mut min = Self::Point::from_value(Self::Scalar::zero());
        let mut max = min;
        for i in 0..Self::Point::len() {
            min[i] = ::cmp_min(l[i], r[i]);
            max[i] = ::cmp_max(l[i], r[i]);
        }

        BoundBox::from_bounds(min, max)
    }
}

impl<P> Linear for Ray<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom + BaseFloat,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise + ElementWise<P!(::Scalar)>
{
    #[inline]
    fn bounding_box(&self) -> BoundBox<P> {
        let outer_bound = Self::Point::from_vec(self.dir()).mul_element_wise(Self::Scalar::infinity());

        let (l, r) = (self.origin, outer_bound);
        let mut min = Self::Point::from_value(Self::Scalar::zero());
        let mut max = min;
        for i in 0..Self::Point::len() {
            min[i] = ::cmp_min(l[i], r[i]);
            max[i] = ::cmp_max(l[i], r[i]);
        }

        BoundBox::from_bounds(min, max)
    }
}

impl<P> Linear for Segment<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise
{
    type Scalar = P::Scalar;
    type Point = P;
    type Vector = P::Diff;

    #[inline]
    fn origin(&self) -> P {
        self.start
    }
    #[inline]
    fn dir(&self) -> P::Diff {
        self.end - self.start
    }
    #[inline]
    fn start(&self) -> Option<P> {
        Some(self.start)
    }
    #[inline]
    fn end(&self) -> Option<P> {
        Some(self.end)
    }
}

impl<P> Linear for Line<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise
{
    type Scalar = P::Scalar;
    type Point = P;
    type Vector = P::Diff;

    fn origin(&self) -> P {
        self.origin
    }

    fn dir(&self) -> P::Diff {
        self.dir
    }

    fn start(&self) -> Option<P> {None}
    fn end(&self) -> Option<P> {None}
}

impl<P: EuclideanSpace> Segment<P> {
    #[inline]
    pub fn new(start: P, end: P) -> Segment<P> {
        Segment{ start, end }
    }
}

macro_rules! inherent_impl_segment {
    ($PointN:ident, $VectorN:ident; $new:ident; ($($start:ident),+), ($($end:ident),+)) => {
        impl<S: BaseScalarGeom> Segment<$PointN<S>> {
            #[inline]
            pub fn $new($($start: S),+, $($end: S),+) -> Segment<$PointN<S>> {
                Segment {
                    start: $PointN::new($($start),+),
                    end: $PointN::new($($end),+)
                }
            }
        }
    }
}
macro_rules! inherent_impl_ray {
    ($Name:ident; $PointN:ident, $VectorN:ident; $new:ident; ($($origin:ident),+), ($($dir:ident),+)) => {
        impl<S: BaseScalarGeom> $Name<$PointN<S>> {
            #[inline]
            pub fn $new($($origin: S),+, $($dir: S),+) -> $Name<$PointN<S>> {
                $Name {
                    origin: $PointN::new($($origin),+),
                    dir: $VectorN::new($($dir),+)
                }
            }
        }
    }
}

inherent_impl_segment!(Point1, Vector1; new1; (start_x), (end_x));
inherent_impl_segment!(Point2, Vector2; new2; (start_x, start_y), (end_x, end_y));
inherent_impl_segment!(Point3, Vector3; new3; (start_x, start_y, start_z), (end_x, end_y, end_z));
inherent_impl_ray!(Ray; Point1, Vector1; new1; (origin_x), (dir_x));
inherent_impl_ray!(Ray; Point2, Vector2; new2; (origin_x, origin_y), (dir_x, dir_y));
inherent_impl_ray!(Ray; Point3, Vector3; new3; (origin_x, origin_y, origin_z), (dir_x, dir_y, dir_z));
inherent_impl_ray!(Line; Point1, Vector1; new1; (origin_x), (dir_x));
inherent_impl_ray!(Line; Point2, Vector2; new2; (origin_x, origin_y), (dir_x, dir_y));
inherent_impl_ray!(Line; Point3, Vector3; new3; (origin_x, origin_y, origin_z), (dir_x, dir_y, dir_z));

impl<L, R> Intersect<R> for L
    where L::Scalar: BaseScalarGeom,
          R: Linear<Point=L::Point, Scalar=L::Scalar, Vector=L::Vector>,
          L: Linear<Point=Point2<<L as Linear>::Scalar>, Vector=Vector2<<L as Linear>::Scalar>>
{
    type Intersection = Point2<L::Scalar>;
    default fn intersect(self, rhs: R) -> Intersection<Point2<L::Scalar>> {
        let (lo, ro) = (self.origin(), rhs.origin());
        let (ld, rd) = (self.dir(), rhs.dir());
        let slope_diff = rd.y*ld.x - rd.x*ld.y;

        if slope_diff == L::Scalar::zero() {
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

impl<P: EuclideanSpace> PartialEq for Line<P>
    where P::Diff: PartialEq + Mul + Array<Element=P::Scalar>,
          P::Scalar: Mul
{
    #[inline]
    default fn eq(&self, other: &Line<P>) -> bool {
        self.dir == other.dir && (other.origin - self.origin).product() == self.dir.product()
    }
}

impl<P: EuclideanSpace> PartialEq for Line<P>
    where P::Diff: PartialEq + Mul + Array<Element=P::Scalar>,
          P::Scalar: Mul + Signed
{
    #[inline]
    fn eq(&self, other: &Line<P>) -> bool {
        self.dir == other.dir && (other.origin - self.origin).product().abs() == self.dir.product().abs()
    }
}
impl<P> Eq for Line<P>
    where P: EuclideanSpace,
          P::Diff: PartialEq + Mul + Array<Element=P::Scalar>,
          P::Scalar: Mul + Eq {}

