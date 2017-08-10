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
pub struct BoundRectTLO<S> {
    pub tl: Point2<S>,
    pub br: Point2<S>
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BoundRectBLO<S> {
    pub bl: Point2<S>,
    pub tr: Point2<S>
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

pub trait RectangleTLO: Rectangle {
    #[inline]
    fn top_left_tlo(&self) -> Point2<Self::Scalar> {
        self.origin_corner()
    }
    #[inline]
    fn top_right_tlo(&self) -> Point2<Self::Scalar> {
        let mut tr = self.origin_corner();
        tr.x += self.width();
        tr
    }
    #[inline]
    fn bottom_left_tlo(&self) -> Point2<Self::Scalar> {
        let mut bl = self.origin_corner();
        bl.y += self.height();
        bl
    }
    #[inline]
    fn bottom_right_tlo(&self) -> Point2<Self::Scalar> {
        self.origin_corner() + self.dims()
    }
}

pub trait RectangleBLO: Rectangle {
    #[inline]
    fn top_left_blo(&self) -> Point2<Self::Scalar> {
        let mut tl = self.origin_corner();
        tl.y += self.height();
        tl
    }
    #[inline]
    fn top_right_blo(&self) -> Point2<Self::Scalar> {
        self.origin_corner() + self.dims()
    }
    #[inline]
    fn bottom_left_blo(&self) -> Point2<Self::Scalar> {
        self.origin_corner()
    }
    #[inline]
    fn bottom_right_blo(&self) -> Point2<Self::Scalar> {
        let mut br = self.origin_corner();
        br.x += self.width();
        br
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

impl<S> BoundRectTLO<S> {
    #[inline]
    pub fn new(top: S, left: S, bottom: S, right: S) -> BoundRectTLO<S> {
        BoundRectTLO {
            tl: Point2{ x: top, y: left },
            br: Point2{ x: bottom, y: right }
        }
    }

    #[inline]
    pub fn cast<T>(&self) -> Option<BoundRectTLO<T>>
        where T: NumCast,
              S: ToPrimitive + Copy
    {
        T::from(self.tl.x)
            .and_then(|x| T::from(self.tl.y).map(|y| Point2{ x, y }))
            .and_then(|tl| T::from(self.br.x).map(|x| (tl, x)))
            .and_then(|(tl, x)| T::from(self.br.y).map(|y| (tl, Point2{ x, y })))
            .map(|(tl, br)| BoundRectTLO{ tl, br })
    }
}

impl<S> BoundRectBLO<S> {
    #[inline]
    pub fn new(bottom: S, left: S, top: S, right: S) -> BoundRectBLO<S> {
        BoundRectBLO {
            bl: Point2{ x: bottom, y: left },
            tr: Point2{ x: top, y: right }
        }
    }

    #[inline]
    pub fn cast<T>(&self) -> Option<BoundRectBLO<T>>
        where T: NumCast,
              S: ToPrimitive + Copy
    {
        T::from(self.bl.x)
            .and_then(|x| T::from(self.bl.y).map(|y| Point2{ x, y }))
            .and_then(|bl| T::from(self.tr.x).map(|x| (bl, x)))
            .and_then(|(bl, x)| T::from(self.tr.y).map(|y| (bl, Point2{ x, y })))
            .map(|(bl, tr)| BoundRectBLO{ bl, tr })
    }
}

impl<S: BaseNum> RectangleTLO for DimsRect<S> {}
impl<S: BaseNum> RectangleBLO for DimsRect<S> {}
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

impl<S: BaseNum> RectangleTLO for OffsetRect<S> {}
impl<S: BaseNum> RectangleBLO for OffsetRect<S> {}
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

impl<S: BaseNum> RectangleTLO for BoundRectTLO<S> {}
impl<S: BaseNum> Rectangle for BoundRectTLO<S> {
    type Scalar = S;

    #[inline]
    fn origin_corner(&self) -> Point2<S> {
        self.tl
    }
    #[inline]
    fn dims(&self) -> Vector2<S> {
        self.br.to_vec() - self.tl.to_vec()
    }
}

impl<S: BaseNum> RectangleBLO for BoundRectBLO<S> {}
impl<S: BaseNum> Rectangle for BoundRectBLO<S> {
    type Scalar = S;

    #[inline]
    fn origin_corner(&self) -> Point2<S> {
        self.bl
    }
    #[inline]
    fn dims(&self) -> Vector2<S> {
        self.tr.to_vec() - self.bl.to_vec()
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

impl<S: BaseNum> From<DimsRect<S>> for BoundRectBLO<S> {
    #[inline]
    fn from(rect: DimsRect<S>) -> BoundRectBLO<S> {
        BoundRectBLO {
            bl: Point2::from_value(S::zero()),
            tr: Point2::from_vec(rect.dims)
        }
    }
}

impl<S: BaseNum> From<DimsRect<S>> for BoundRectTLO<S> {
    #[inline]
    fn from(rect: DimsRect<S>) -> BoundRectTLO<S> {
        BoundRectTLO {
            tl: Point2::from_value(S::zero()),
            br: Point2::from_vec(rect.dims)
        }
    }
}

impl<S: BaseNum> From<OffsetRect<S>> for BoundRectBLO<S> {
    #[inline]
    fn from(rect: OffsetRect<S>) -> BoundRectBLO<S> {
        BoundRectBLO {
            bl: rect.origin,
            tr: rect.origin + rect.dims
        }
    }
}

impl<S: BaseNum> From<OffsetRect<S>> for BoundRectTLO<S> {
    #[inline]
    fn from(rect: OffsetRect<S>) -> BoundRectTLO<S> {
        BoundRectTLO {
            tl: rect.origin,
            br: rect.origin + rect.dims
        }
    }
}

impl<S: BaseNum> From<BoundRectBLO<S>> for OffsetRect<S> {
    #[inline]
    fn from(rect: BoundRectBLO<S>) -> OffsetRect<S> {
        OffsetRect {
            origin: rect.bl,
            dims: rect.tr.to_vec() - rect.bl.to_vec()
        }
    }
}

impl<S: BaseNum> From<BoundRectTLO<S>> for OffsetRect<S> {
    #[inline]
    fn from(rect: BoundRectTLO<S>) -> OffsetRect<S> {
        OffsetRect {
            origin: rect.tl,
            dims: rect.br.to_vec() - rect.tl.to_vec()
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

impl<S: BaseNum> Add<Vector2<S>> for BoundRectBLO<S> {
    type Output = Self;
    #[inline]
    fn add(mut self, rhs: Vector2<S>) -> BoundRectBLO<S> {
        self.bl += rhs;
        self.tr += rhs;
        self
    }
}

impl<S: BaseNum> Add<Vector2<S>> for BoundRectTLO<S> {
    type Output = Self;
    #[inline]
    fn add(mut self, rhs: Vector2<S>) -> BoundRectTLO<S> {
        self.tl += rhs;
        self.br += rhs;
        self
    }
}

impl<S: BaseNum> Sub<Vector2<S>> for BoundRectBLO<S> {
    type Output = Self;
    #[inline]
    fn sub(mut self, rhs: Vector2<S>) -> BoundRectBLO<S> {
        self.bl -= rhs;
        self.tr -= rhs;
        self
    }
}

impl<S: BaseNum> Sub<Vector2<S>> for BoundRectTLO<S> {
    type Output = Self;
    #[inline]
    fn sub(mut self, rhs: Vector2<S>) -> BoundRectTLO<S> {
        self.tl -= rhs;
        self.br -= rhs;
        self
    }
}
