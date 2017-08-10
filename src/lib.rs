extern crate cgmath;
extern crate num_traits;

use cgmath::*;
use std::ops::{Add, Sub};
use num_traits::{Bounded, NumCast, ToPrimitive};

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
    type Scalar: BaseNum;

    fn origin_corner(&self) -> Point2<Self::Scalar>;
    fn dims(&self) -> Vector2<Self::Scalar>;

    #[inline]
    fn contains(&self, point: Point2<Self::Scalar>) -> bool {
        let origin_corner = self.origin_corner();
        let dims = self.dims();

        origin_corner.x <= point.x &&
        origin_corner.y <= point.y &&
        point.x <= origin_corner.x + dims.x &&
        point.y <= origin_corner.y + dims.y
    }

    #[inline]
    fn width(&self) -> Self::Scalar {
        self.dims().x
    }
    #[inline]
    fn height(&self) -> Self::Scalar {
        self.dims().y
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

impl<S: BaseNum> Rectangle for DimsRect<S> {
    type Scalar = S;

    #[inline]
    fn origin_corner(&self) -> Point2<S> { Point2::new(S::zero(), S::zero())
    }
    #[inline]
    fn dims(&self) -> Vector2<S> {
        self.dims
    }
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

impl<S: BaseNum> Rectangle for OffsetRect<S> {
    type Scalar = S;

    #[inline]
    fn origin_corner(&self) -> Point2<S> {
        self.origin
    }
    #[inline]
    fn dims(&self) -> Vector2<S> {
        self.dims
    }
}

impl<S: BaseNum> Rectangle for BoundRect<S> {
    type Scalar = S;

    #[inline]
    fn origin_corner(&self) -> Point2<S> {
        self.min
    }
    #[inline]
    fn dims(&self) -> Vector2<S> {
        self.max.to_vec() - self.min.to_vec()
    }
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
