use {MulDiv, BaseNumGeom};
use cgmath::*;

use std::fmt::Debug;

use line::Line;

use std::ops::{Add, Sub};
use num_traits::{Bounded, NumCast, ToPrimitive, PrimInt};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DimsRect<S> {
    pub dims: Vector2<S>
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct OffsetRect<S> {
    pub origin: Point2<S>,
    pub dims: Vector2<S>
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BoundRect<S> {
    pub min: Point2<S>,
    pub max: Point2<S>
}

pub trait Rectangle {
    type Scalar: BaseNum + MulDiv;
    type Point: EuclideanSpace<Scalar=Self::Scalar, Diff=Self::Vector> + ElementWise<Self::Scalar> + MulDiv<Self::Scalar> + Debug;
    type Vector: VectorSpace<Scalar=Self::Scalar> + Array<Element=Self::Scalar> + MulDiv + MulDiv<Self::Scalar> + PartialEq + Debug;

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
        self.dims()[0]
    }
    #[inline]
    fn height(&self) -> Self::Scalar {
        self.dims()[1]
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

    #[inline]
    fn intersects_int<L>(&self, line: L) -> (Option<Self::Point>, Option<Self::Point>)
        where L: Line<Scalar=Self::Scalar, Point=Self::Point, Vector=Self::Vector>,
              Self::Scalar: PrimInt + Debug
    {
        let min = self.min();
        let max = self.max();
        let zero = Self::Scalar::zero();

        let (start, end) = (line.start(), line.end());
        let (mut enter, mut exit) = (start, end);
        let (mut enter_valid, mut exit_valid) = (false, false);
        let dir = line.dir();

        for i in 0..Self::Point::len() {
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

        for i in 0..Self::Point::len() {
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

impl<S> DimsRect<S> {
    #[inline]
    pub fn new(width: S, height: S) -> DimsRect<S> {
        DimsRect {
            dims: Vector2::new(width, height)
        }
    }

    #[inline]
    pub fn cast<T>(&self) -> Option<DimsRect<T>>
        where T: NumCast,
              S: ToPrimitive + Copy
    {
        T::from(self.dims.x)
            .and_then(|x| T::from(self.dims.y).map(|y| Vector2::new(x, y)))
            .map(|dims| DimsRect{ dims })
    }
}

impl<S> OffsetRect<S> {
    #[inline]
    pub fn new(origin_x: S, origin_y: S, width: S, height: S) -> OffsetRect<S> {
        OffsetRect {
            origin: Point2{ x: origin_x, y: origin_y },
            dims: Vector2::new(width, height)
        }
    }

    #[inline]
    pub fn cast<T>(&self) -> Option<OffsetRect<T>>
        where T: NumCast,
              S: ToPrimitive + Copy
    {
        T::from(self.origin.x)
            .and_then(|x| T::from(self.origin.y).map(|y| Point2{ x, y }))
            .and_then(|origin| T::from(self.dims.x).map(|x| (origin, x)))
            .and_then(|(origin, x)| T::from(self.dims.y).map(|y| (origin, Vector2{ x, y })))
            .map(|(origin, dims)| OffsetRect{ origin, dims })
    }
}

impl<S> BoundRect<S> {
    #[inline]
    pub fn new(min_x: S, min_y: S, max_x: S, max_y: S) -> BoundRect<S> {
        BoundRect {
            min: Point2{ x: min_x, y: min_y },
            max: Point2{ x: max_x, y: max_y }
        }
    }

    #[inline]
    pub fn cast<T>(&self) -> Option<BoundRect<T>>
        where T: NumCast,
              S: ToPrimitive + Copy
    {
        T::from(self.min.x)
            .and_then(|x| T::from(self.min.y).map(|y| Point2{ x, y }))
            .and_then(|min| T::from(self.max.x).map(|x| (min, x)))
            .and_then(|(min, x)| T::from(self.max.y).map(|y| (min, Point2{ x, y })))
            .map(|(min, max)| BoundRect{ min, max })
    }
}

impl<S: BaseNumGeom> Rectangle for DimsRect<S> {
    type Scalar = S;
    type Point = Point2<S>;
    type Vector = Vector2<S>;

    #[inline]
    fn min(&self) -> Point2<S> {Point2::new(S::zero(), S::zero())}
    #[inline]
    fn dims(&self) -> Vector2<S> {self.dims}
}

impl<S: Bounded> Bounded for DimsRect<S> {
    #[inline]
    fn min_value() -> DimsRect<S> {
        DimsRect {
            dims: Vector2::min_value()
        }
    }

    #[inline]
    fn max_value() -> DimsRect<S> {
        DimsRect {
            dims: Vector2::max_value()
        }
    }
}

impl<S: BaseNumGeom> Rectangle for OffsetRect<S> {
    type Scalar = S;
    type Point = Point2<S>;
    type Vector = Vector2<S>;

    #[inline]
    fn min(&self) -> Point2<S> {self.origin}
    #[inline]
    fn dims(&self) -> Vector2<S> {self.dims}
}

impl<S: BaseNumGeom> Rectangle for BoundRect<S> {
    type Scalar = S;
    type Point = Point2<S>;
    type Vector = Vector2<S>;

    #[inline]
    fn min(&self) -> Point2<S> {self.min}
    #[inline]
    fn max(&self) -> Point2<S> {self.max}
}

impl<S: BaseNum> From<DimsRect<S>> for OffsetRect<S> {
    #[inline]
    fn from(rect: DimsRect<S>) -> OffsetRect<S> {
        OffsetRect {
            origin: Point2::from_value(S::zero()),
            dims: rect.dims
        }
    }
}

impl<S: BaseNum> From<DimsRect<S>> for BoundRect<S> {
    #[inline]
    fn from(rect: DimsRect<S>) -> BoundRect<S> {
        BoundRect {
            min: Point2::from_value(S::zero()),
            max: Point2::from_vec(rect.dims)
        }
    }
}

impl<S: BaseNum> From<OffsetRect<S>> for BoundRect<S> {
    #[inline]
    fn from(rect: OffsetRect<S>) -> BoundRect<S> {
        BoundRect {
            min: rect.origin,
            max: rect.origin + rect.dims
        }
    }
}

impl<S: BaseNum> From<BoundRect<S>> for OffsetRect<S> {
    #[inline]
    fn from(rect: BoundRect<S>) -> OffsetRect<S> {
        OffsetRect {
            origin: rect.min,
            dims: rect.max.to_vec() - rect.min.to_vec()
        }
    }
}

impl<S: BaseNum> Add<Vector2<S>> for OffsetRect<S> {
    type Output = Self;
    #[inline]
    fn add(mut self, rhs: Vector2<S>) -> OffsetRect<S> {
        self.origin += rhs;
        self
    }
}

impl<S: BaseNum> Sub<Vector2<S>> for OffsetRect<S> {
    type Output = Self;
    #[inline]
    fn sub(mut self, rhs: Vector2<S>) -> OffsetRect<S> {
        self.origin -= rhs;
        self
    }
}

impl<S: BaseNum> Add<Vector2<S>> for BoundRect<S> {
    type Output = Self;
    #[inline]
    fn add(mut self, rhs: Vector2<S>) -> BoundRect<S> {
        self.min += rhs;
        self.max += rhs;
        self
    }
}

impl<S: BaseNum> Sub<Vector2<S>> for BoundRect<S> {
    type Output = Self;
    #[inline]
    fn sub(mut self, rhs: Vector2<S>) -> BoundRect<S> {
        self.min -= rhs;
        self.max -= rhs;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use line::Segment;

    #[test]
    fn line_touch_edge_exit() {
        let rect = BoundRect::new(20, 20, 30, 40);
        let segment = Segment::new(25, 25, 20, 25);
        assert_eq!((None, None), rect.intersects_int(segment));
    }

    #[test]
    fn line_touch_edge_outer_exit() {
        let rect = BoundRect::new(20, 20, 30, 40);
        let segment = Segment::new(25, 15, 25, 20);
        assert_eq!((Some(Point2::new(25, 20)), None), rect.intersects_int(segment));
    }

    #[test]
    fn line_exit_single_move_lt() {
        let rect = BoundRect::new(20, 20, 30, 40);
        let segment = Segment::new(20, 25, 19, 25);
        assert_eq!((None, Some(Point2::new(20, 25))), rect.intersects_int(segment));
    }

    #[test]
    fn line_exit_single_move_gt() {
        let rect = BoundRect::new(20, 20, 30, 40);
        let segment = Segment::new(30, 25, 31, 25);
        assert_eq!((None, Some(Point2::new(30, 25))), rect.intersects_int(segment));
    }

    #[test]
    fn line_touch_edge_enter() {
        let rect = BoundRect::new(20, 20, 30, 40);
        let segment = Segment::new(25, 0, 27, 20);
        assert_eq!((Some(Point2::new(27, 20)), None), rect.intersects_int(segment));
    }

    #[test]
    fn intersect_diagonal() {
        let rect = BoundRect::new(20, 20, 40, 40);
        let segment = Segment::new(0, 0, 50, 50);
        assert_eq!((Some(Point2::new(20, 20)), Some(Point2::new(40, 40))), rect.intersects_int(segment));
    }

    #[test]
    fn intersect_horizontal() {
        let rect = BoundRect::new(20, 20, 30, 40);
        let segment = Segment::new(0, 25, 50, 25);
        assert_eq!((Some(Point2::new(20, 25)), Some(Point2::new(30, 25))), rect.intersects_int(segment));
    }
}
