extern crate cgmath;
extern crate num_traits;

use cgmath::*;
use num_traits::Bounded;

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

impl<S: BaseNum> RectangleTLO for DimsRect<S> {}
impl<S: BaseNum> RectangleBLO for DimsRect<S> {}
impl<S: BaseNum> Rectangle for DimsRect<S> {
    type Scalar = S;

    #[inline]
    fn origin_corner(&self) -> Point2<S> {
        Point2::new(S::zero(), S::zero())
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
