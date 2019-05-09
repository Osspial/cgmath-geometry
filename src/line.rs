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
use num_traits::{Bounded, Float};
use std::cmp::Ordering;
use rect::{BoundBox, GeoBox};

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct Ray<D: Dimensionality<S>, S: BaseScalarGeom> {
    pub origin: D::Point,
    pub dir: D::Vector
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct Segment<D: Dimensionality<S>, S: BaseScalarGeom> {
    pub start: D::Point,
    pub end: D::Point
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature="serde", derive(Deserialize, Serialize))]
pub struct Line<D: Dimensionality<S>, S: BaseScalarGeom> {
    pub origin: D::Point,
    pub dir: D::Vector
}

impl<D: Dimensionality<S>, S: BaseScalarGeom> Segment<D, S> {
    pub fn bounding_box(&self) -> BoundBox<D, S> {
        let l = self.start;
        let r = self.end;

        let mut min = D::Point::from_value(S::zero());
        let mut max = min;
        for i in 0..D::Point::len() {
            min[i] = ::cmp_min(l[i], r[i]);
            max[i] = ::cmp_max(l[i], r[i]);
        }

        BoundBox::from_bounds(min, max)
    }
}

pub trait Linear
    where d!(Vector): VectorSpace<Scalar=Self::Scalar>,
          d!(Point): EuclideanSpace<Scalar=Self::Scalar, Diff=d!(Vector)>
{
    type Scalar: BaseScalarGeom;
    type D: Dimensionality<Self::Scalar>;

    fn origin(&self) -> d!(Point);
    fn dir(&self) -> d!(Vector);
    #[inline]
    fn dir_recip(&self) -> d!(Vector)
        where Self::Scalar: Float
    {
        let mut d = self.dir();
        for i in 0..d!(Vector::len()) {
            d[i] = d[i].recip();
        }
        d
    }
    fn start(&self) -> Option<d!(Point)>;
    fn end(&self) -> Option<d!(Point)>;

    fn clip_to_scalar_bounds(&self) -> Segment<Self::D, Self::Scalar> {
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

            for i in 0..d!(Point::len()) {
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

impl<D, S> Linear for Ray<D, S>
    where S: BaseScalarGeom,
          D: Dimensionality<S>,
          D::Vector: VectorSpace<Scalar=S>,
          D::Point: EuclideanSpace<Scalar=S, Diff=D::Vector>
{
    type Scalar = S;
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
}

impl<D, S> Linear for Segment<D, S>
    where S: BaseScalarGeom,
          D: Dimensionality<S>,
          D::Vector: VectorSpace<Scalar=S>,
          D::Point: EuclideanSpace<Scalar=S, Diff=D::Vector>
{
    type Scalar = S;
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

impl<D, S> Linear for Line<D, S>
    where S: BaseScalarGeom,
          D: Dimensionality<S>,
          D::Vector: VectorSpace<Scalar=S>,
          D::Point: EuclideanSpace<Scalar=S, Diff=D::Vector>
{
    type Scalar = S;
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

impl<D: Dimensionality<S>, S: BaseScalarGeom> Segment<D, S> {
    #[inline]
    pub fn new(start: D::Point, end: D::Point) -> Segment<D, S> {
        Segment{ start, end }
    }
}

macro_rules! inherent_impl_segment {
    ($D:ident; $new:ident; ($($start:ident),+), ($($end:ident),+)) => {
        impl<S: BaseScalarGeom> Segment<$D, S> {
            #[inline]
            pub fn $new($($start: S),+, $($end: S),+) -> Segment<$D, S> {
                Segment {
                    start: <$D as Dimensionality<S>>::Point::new($($start),+),
                    end: <$D as Dimensionality<S>>::Point::new($($end),+)
                }
            }
        }
    }
}
macro_rules! inherent_impl_ray {
    ($Name:ident; $D:ident; $new:ident; ($($origin:ident),+), ($($dir:ident),+)) => {
        impl<S: BaseScalarGeom> $Name<$D, S> {
            #[inline]
            pub fn $new($($origin: S),+, $($dir: S),+) -> $Name<$D, S> {
                $Name {
                    origin: <$D as Dimensionality<S>>::Point::new($($origin),+),
                    dir: <$D as Dimensionality<S>>::Vector::new($($dir),+)
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
    ($($t:tt)*) => {<L::D as Dimensionality<L::Scalar>>::$($t)* };
}

impl<L, R> Intersect<R> for L
    where L: Linear<D=D2>,
          R: Linear<Scalar=L::Scalar, D=D2>,
{
    type Intersection = Point2<L::Scalar>;
    fn intersect(self, rhs: R) -> Intersection<Point2<L::Scalar>> {
        // To see how this works, check out https://www.desmos.com/calculator/oipyx0c0me

        let (lo, ro) = (self.origin(), rhs.origin());
        let (ld, rd) = (self.dir(), rhs.dir());

        let slope_diff = ld.x * rd.y - ld.y * rd.x;

        if slope_diff == L::Scalar::zero() {
            return if (ro.y-lo.y) * (ro.x-lo.x) == ld.x * ld.y || ro == lo {
                Intersection::Eq
            } else {
                Intersection::None
            };
        }
        // This could be optimized to computing `t` then lerping between lo and (lo + ld) for
        // floats, but that approach breaks down for integers.
        let intersection = Point2::new(
            (ro.x * ld.x * rd.y - rd.x * (lo.x * ld.y - ld.x * (lo.y - ro.y))) / slope_diff,
            (rd.y * (lo.y * ld.x - ld.y * (lo.x - ro.x)) - ro.y * ld.y * rd.x) / slope_diff,
        );

        match bounding_box(self).intersect_rect(bounding_box(rhs)).map(|r| r.contains(intersection)) {
            Some(true) => Intersection::Some(intersection),
            _ => Intersection::None
        }
    }
}

fn bounding_box<L: Linear>(line: L) -> BoundBox<L::D, L::Scalar>
    where L::Scalar: BaseScalarGeom,
          L::D: Dimensionality<L::Scalar>,
          ld!(Vector): VectorSpace<Scalar=L::Scalar>,
          ld!(Point): EuclideanSpace<Scalar=L::Scalar, Diff=ld!(Vector)>
{
    let origin = line.origin();
    let mut end_bound = ld!(Point::from_vec(line.dir()));
    let mut start_bound = ld!(Point::from_vec(line.dir()));
    for i in 0..ld!(Point::len()) {
        match end_bound[i].partial_cmp(&L::Scalar::zero()) {
            Some(Ordering::Less) => {
                end_bound[i] = L::Scalar::min_value();
                start_bound[i] = L::Scalar::max_value();
            },
            Some(Ordering::Greater) => {
                end_bound[i] = L::Scalar::max_value();
                start_bound[i] = L::Scalar::min_value();
            },
            None |
            Some(Ordering::Equal) => {
                end_bound[i] = origin[i];
                start_bound[i] = origin[i];
            },
        }
    }

    let start = line.start().unwrap_or(start_bound);
    let end = line.end().unwrap_or(end_bound);

    let mut min = ld!(Point::from_value(L::Scalar::zero()));
    let mut max = min;
    for i in 0..ld!(Point::len()) {
        min[i] = ::cmp_min(start[i], end[i]);
        max[i] = ::cmp_max(start[i], end[i]);
    }

    BoundBox::from_bounds(min, max)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn line_intersect() {
        let horizontal = Segment::new2(-1.0, 0.0, 1.0, 0.0);
        let vertical = Segment::new2(0.0, -1.0, 0.0, 1.0);
        assert_eq!(Intersection::Some(Point2::new(0.0, 0.0)), horizontal.intersect(vertical));

        let horizontal = Segment::new2(-1.0, 0.0, 1.0, 0.0);
        let vertical = Segment::new2(0.5, -1.0, 0.5, 1.0);
        assert_eq!(Intersection::Some(Point2::new(0.5, 0.0)), horizontal.intersect(vertical));

        let diag0 = Segment::new2(-1.0, -1.0, 1.0, 1.0);
        let diag1 = Segment::new2(-1.0, 1.0, 1.0, -1.0);
        assert_eq!(Intersection::Some(Point2::new(0.0, 0.0)), diag0.intersect(diag1));


        let horizontal = Segment::new2(-1, 0, 1, 0);
        let vertical = Segment::new2(0, -1, 0, 1);
        assert_eq!(Intersection::Some(Point2::new(0, 0)), horizontal.intersect(vertical));

        let diag0 = Segment::new2(-1, -1, 1, 1);
        let diag1 = Segment::new2(-1, 1, 1, -1);
        assert_eq!(Intersection::Some(Point2::new(0, 0)), diag0.intersect(diag1));

        let diag0 = Segment::new2(-1, -1, 2, 1);
        let diag1 = Segment::new2(-1, 1, 2, -1);
        assert_eq!(Intersection::Some(Point2::new(0, 0)), diag0.intersect(diag1));
    }

    #[test]
    fn general_bounding_box() {
        let imin = i32::min_value();
        let imax = i32::max_value();

        let ray: Ray<D2, i32> = Ray {
            origin: Point2::new(50, 50),
            dir: Vector2::new(1, 0)
        };
        let ray_bounds = BoundBox::new2(
            50, 50,
            imax, 50
        );
        assert_eq!(ray_bounds, bounding_box(ray));

        let ray: Ray<D2, i32> = Ray {
            origin: Point2::new(50, 50),
            dir: Vector2::new(1, 2)
        };
        let ray_bounds = BoundBox::new2(
            50, 50,
            imax, imax
        );
        assert_eq!(ray_bounds, bounding_box(ray));

        let ray: Ray<D2, i32> = Ray {
            origin: Point2::new(50, 50),
            dir: Vector2::new(-1, 2)
        };
        let ray_bounds = BoundBox::new2(
            imin, 50,
            50, imax
        );
        assert_eq!(ray_bounds, bounding_box(ray));

        let ray: Ray<D2, i32> = Ray {
            origin: Point2::new(50, 50),
            dir: Vector2::new(-1, -1)
        };
        let ray_bounds = BoundBox::new2(
            imin, imin,
            50, 50
        );
        assert_eq!(ray_bounds, bounding_box(ray));

        let ray: Ray<D2, i32> = Ray {
            origin: Point2::new(50, 50),
            dir: Vector2::new(0, 0)
        };
        let ray_bounds = BoundBox::new2(
            50, 50,
            50, 50
        );
        assert_eq!(ray_bounds, bounding_box(ray));

        // TODO: TEST LINE AND SEGMENT
    }
}
