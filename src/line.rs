use {MulDiv, BaseScalarGeom, Intersect, Intersection};
use cgmath::*;
use num_traits::{Bounded, Float};
use std::cmp::Ordering;
use rect::{BoundBox, Box};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Ray<P: EuclideanSpace> {
    pub origin: P,
    pub dir: P::Diff
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Segment<P: EuclideanSpace> {
    pub start: P,
    pub end: P
}

pub trait Line {
    type Scalar: BaseScalarGeom;
    type Point: EuclideanSpace<Scalar=Self::Scalar, Diff=Self::Vector> + ElementWise<Self::Scalar> + MulDiv<Self::Scalar>;
    type Vector: VectorSpace<Scalar=Self::Scalar> + Array<Element=Self::Scalar> + MulDiv + MulDiv<Self::Scalar> + ElementWise;

    fn origin(&self) -> Self::Point;
    fn dir(&self) -> Self::Vector;
    fn start(&self) -> Self::Point;
    fn end(&self) -> Self::Point;

    fn bounding_box(&self) -> BoundBox<Self::Point> {
        let l = self.start();
        let r = self.end();

        let mut min = Self::Point::from_value(Self::Scalar::zero());
        let mut max = min;
        for i in 0..Self::Point::len() {
            min[i] = ::cmp_min(l[i], r[i]);
            max[i] = ::cmp_max(l[i], r[i]);
        }

        BoundBox::from_bounds(min, max)
    }
}

impl<P> Line for Ray<P>
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
    fn start(&self) -> P {
        self.origin
    }
    #[inline]
    default fn end(&self) -> P {
        let origin = self.origin();
        let dir = self.dir();
        let zero = Self::Scalar::zero();

        let mut axis_closest = usize::max_value();
        let mut axis_closest_threshold = Self::Scalar::max_value();
        let mut axis_closest_muldiv = (zero, zero);

        for i in 0..Self::Point::len() {
            let axis_dist = match dir[i].partial_cmp(&Self::Scalar::zero()) {
                Some(Ordering::Less) => zero - (origin[i] - Self::Scalar::min_value()),
                Some(Ordering::Greater) => Self::Scalar::max_value() - origin[i],
                Some(Ordering::Equal) |
                None                 => Self::Scalar::max_value()
            };

            let axis_threshold = axis_dist / dir[i];
            if axis_threshold <= axis_closest_threshold {
                axis_closest = i;
                axis_closest_threshold = axis_threshold;
                axis_closest_muldiv = (axis_dist, dir[i]);
            }
        }

        if axis_closest == usize::max_value() {
            panic!("direction equal to zero");
        } else {
            origin + dir.mul_div(axis_closest_muldiv.0, axis_closest_muldiv.1)
        }
    }
}

impl<P> Line for Ray<P>
    where P: EuclideanSpace + ElementWise<P!(::Scalar)> + MulDiv<P!(::Scalar)>,
          P!(::Scalar): BaseScalarGeom + BaseFloat,
          P::Diff: VectorSpace<Scalar=P!(::Scalar)> + Array<Element=P!(::Scalar)> + MulDiv + MulDiv<P!(::Scalar)> + ElementWise + ElementWise<P!(::Scalar)>
{
    #[inline]
    fn end(&self) -> P {
        P::from_vec(self.dir.mul_element_wise(<P::Scalar as Float>::infinity()))
    }
}

impl<P> Line for Segment<P>
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
    fn start(&self) -> P {
        self.start
    }
    #[inline]
    fn end(&self) -> P {
        self.end
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
    ($PointN:ident, $VectorN:ident; $new:ident; ($($origin:ident),+), ($($dir:ident),+)) => {
        impl<S: BaseScalarGeom> Ray<$PointN<S>> {
            #[inline]
            pub fn $new($($origin: S),+, $($dir: S),+) -> Ray<$PointN<S>> {
                Ray {
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
inherent_impl_ray!(Point1, Vector1; new1; (origin_x), (dir_x));
inherent_impl_ray!(Point2, Vector2; new2; (origin_x, origin_y), (dir_x, dir_y));
inherent_impl_ray!(Point3, Vector3; new3; (origin_x, origin_y, origin_z), (dir_x, dir_y, dir_z));

impl<L, R> Intersect<R> for L
    where L::Scalar: BaseScalarGeom,
          R: Line<Point=L::Point, Scalar=L::Scalar, Vector=L::Vector>,
          L: Line<Point=Point2<<L as Line>::Scalar>, Vector=Vector2<<L as Line>::Scalar>>
{
    type Intersection = Point2<L::Scalar>;
    fn intersect(self, rhs: R) -> Intersection<Point2<L::Scalar>> {
        let (lo, ro) = (self.origin(), rhs.origin());
        let (ld, rd) = (self.dir(), rhs.dir());
        let slope_diff = ld.x*rd.y + rd.x*ld.y;

        if slope_diff == L::Scalar::zero() {
            return if (ro.y-lo.y) * (ro.x-lo.x) == ld.x * ld.y || ro == lo {
                Intersection::Eq
            } else {
                Intersection::None
            };
        }
        let t = (rd.y*(lo.x-ro.x) + rd.x*(ro.y-lo.y))/slope_diff;
        let intersection = lo + ld * t;

        match self.bounding_box().intersect_rect(rhs.bounding_box()).map(|r| r.contains(intersection)) {
            Some(true) => Intersection::Some(intersection),
            _ => Intersection::None
        }
    }
}
