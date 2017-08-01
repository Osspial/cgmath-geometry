extern crate cgmath;
extern crate num_traits;

use cgmath::*;
use num_traits::*;
use std::ops::Sub;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Rect<S> {
    pub dims: Vector2<S>
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct OffsetRect<S> {
    pub origin: Point2<S>,
    pub dims: Vector2<S>
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BoundRect<S> {
    pub tl: Point2<S>,
    pub br: Point2<S>
}

pub trait Rectangle {
    type Scalar: Sized;
    fn top_left(&self) -> Point2<Self::Scalar>;
    fn bottom_right(&self) -> Point2<Self::Scalar>;

    #[inline]
    fn width(&self) -> Self::Scalar
        where Self::Scalar: Sub<Output=Self::Scalar>
    {
        self.bottom_right().x - self.top_left().x
    }

    #[inline]
    fn height(&self) -> Self::Scalar
        where Self::Scalar: Sub<Output=Self::Scalar>
    {
        self.bottom_right().y - self.top_left().y
    }
}
