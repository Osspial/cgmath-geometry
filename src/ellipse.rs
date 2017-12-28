use {BaseNumGeom, MulDiv};
use rect::DimsRect;
use cgmath::*;

pub struct Circle2<S> {
    pub origin: Point2<S>,
    pub radius: S
}

pub struct Ellipse2<S> {
    pub origin: Point2<S>,
    pub dims: DimsRect<S>
}

pub trait Ellipse {
    type Scalar: BaseNum + MulDiv;
    type Point: EuclideanSpace<Scalar=Self::Scalar, Diff=Self::Vector> + ElementWise<Self::Scalar> + MulDiv<Self::Scalar>;
    type Vector: VectorSpace<Scalar=Self::Scalar> + Array<Element=Self::Scalar> + MulDiv + MulDiv<Self::Scalar>;

    fn center(&self) -> Self::Point;
    fn dims(&self) -> Self::Vector;

    fn width(&self) -> Self::Scalar {
        self.dims()[0]
    }
    fn height(&self) -> Self::Scalar {
        self.dims()[1]
    }

    fn contains(&self, point: Self::Point) -> bool {
        let mut sum = Self::Scalar::zero();
        let center = self.center();
        let dims = self.dims();

        for i in 0..Self::Point::len() {
            let n = (point[i] - center[i]) / dims[i];
            sum += n*n;
        }

        sum <= Self::Scalar::one()
    }
}

pub trait Circle: Ellipse {
    fn radius(&self) -> Self::Scalar {
        self.dims()[0]
    }
}

impl<S: BaseNumGeom> Ellipse for Circle2<S> {
    type Scalar = S;
    type Point = Point2<S>;
    type Vector = Vector2<S>;

    #[inline(always)]
    fn center(&self) -> Point2<S> {
        self.origin
    }

    #[inline(always)]
    fn dims(&self) -> Vector2<S> {
        Vector2::new(self.radius, self.radius)
    }
}
impl<S: BaseNumGeom> Circle for Circle2<S> {}

impl<S: BaseNumGeom> Ellipse for Ellipse2<S> {
    type Scalar = S;
    type Point = Point2<S>;
    type Vector = Vector2<S>;

    #[inline(always)]
    fn center(&self) -> Point2<S> {
        self.origin
    }

    #[inline(always)]
    fn dims(&self) -> Vector2<S> {
        self.dims.dims
    }
}
