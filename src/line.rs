use {MulDiv, BaseNumGeom};
use cgmath::*;
use num_traits::Bounded;
use std::cmp::Ordering;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Ray<S> {
    pub origin: Point2<S>,
    pub dir: Vector2<S>
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Segment<S> {
    pub start: Point2<S>,
    pub end: Point2<S>
}

pub trait Line {
    type Scalar: BaseNumGeom;
    type Point: EuclideanSpace<Scalar=Self::Scalar, Diff=Self::Vector>;
    type Vector: VectorSpace<Scalar=Self::Scalar> + Array<Element=Self::Scalar>;

    fn origin(&self) -> Self::Point;
    fn dir(&self) -> Self::Vector;
    fn start(&self) -> Self::Point;
    fn end(&self) -> Self::Point;
}

impl<S: BaseNumGeom> Line for Ray<S> {
    type Scalar = S;
    type Point = Point2<S>;
    type Vector = Vector2<S>;

    #[inline]
    fn origin(&self) -> Point2<S> {self.origin}
    #[inline]
    fn dir(&self) -> Vector2<S> {self.dir}
    #[inline]
    fn start(&self) -> Point2<S> {
        self.origin
    }
    #[inline]
    default fn end(&self) -> Point2<S> {
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

impl<F: BaseNumGeom + BaseFloat> Line for Ray<F> {
    #[inline]
    fn end(&self) -> Point2<F> {
        Point2::from_vec(self.dir.mul_element_wise(F::infinity()))
    }
}

impl<S: BaseNumGeom> Line for Segment<S> {
    type Scalar = S;
    type Point = Point2<S>;
    type Vector = Vector2<S>;

    #[inline]
    fn origin(&self) -> Point2<S> {
        self.start
    }
    #[inline]
    fn dir(&self) -> Vector2<S> {
        self.end - self.start
    }
    #[inline]
    fn start(&self) -> Point2<S> {
        self.start
    }
    #[inline]
    fn end(&self) -> Point2<S> {
        self.end
    }
}

impl<S: BaseNumGeom> Segment<S> {
    #[inline]
    pub fn new(start_x: S, start_y: S, end_x: S, end_y: S) -> Segment<S> {
        Segment{ start: Point2::new(start_x, start_y), end: Point2::new(end_x, end_y) }
    }
}

impl<S: BaseNumGeom> Ray<S> {
    #[inline]
    pub fn new(origin_x: S, origin_y: S, dir_x: S, dir_y: S) -> Ray<S> {
        Ray{ origin: Point2::new(origin_x, origin_y), dir: Vector2::new(dir_x, dir_y) }
    }
}
