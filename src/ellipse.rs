use {BaseScalarGeom, MulDiv};
use rect::{DimsBox, Box};
use polar::{Polar2, PolarSpace};

use cgmath::*;
use num_traits::pow::pow;

pub struct Circle2<S> {
    pub origin: Point2<S>,
    pub radius: S
}

pub struct Ellipse2<S: BaseScalarGeom> {
    pub origin: Point2<S>,
    pub dims: DimsBox<Point2<S>>
}

pub trait Ellipse {
    type Scalar: BaseFloat + MulDiv;
    type Point: EuclideanSpace<Scalar=Self::Scalar, Diff=Self::Vector> + ElementWise<Self::Scalar> + MulDiv<Self::Scalar>;
    type Vector: VectorSpace<Scalar=Self::Scalar> + Array<Element=Self::Scalar> + MulDiv + MulDiv<Self::Scalar>;
    type Polar: PolarSpace<Scalar=Self::Scalar>;

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

    fn project_polar(&self, polar: Self::Polar) -> Self::Polar;
}

pub trait Circle: Ellipse {
    fn radius(&self) -> Self::Scalar {
        self.dims()[0]
    }
}

impl<S> Ellipse for Circle2<S>
    where S: BaseFloat + MulDiv
{
    type Scalar = S;
    type Point = Point2<S>;
    type Vector = Vector2<S>;
    type Polar = Polar2<Rad<S>>;

    #[inline(always)]
    fn center(&self) -> Point2<S> {
        self.origin
    }

    #[inline(always)]
    fn dims(&self) -> Vector2<S> {
        Vector2::new(self.radius, self.radius)
    }

    #[inline(always)]
    fn project_polar(&self, mut polar: Polar2<Rad<S>>) -> Polar2<Rad<S>> {
        polar.dist = self.radius;
        polar
    }
}
impl<S> Circle for Circle2<S>
    where S: BaseScalarGeom + BaseFloat {}

impl<S> Ellipse for Ellipse2<S>
    where S: BaseScalarGeom + BaseFloat
{
    type Scalar = S;
    type Point = Point2<S>;
    type Vector = Vector2<S>;
    type Polar = Polar2<Rad<S>>;

    #[inline(always)]
    fn center(&self) -> Point2<S> {
        self.origin
    }

    #[inline(always)]
    fn dims(&self) -> Vector2<S> {
        self.dims.dims
    }

    #[inline(always)]
    fn project_polar(&self, mut polar: Polar2<Rad<S>>) -> Polar2<Rad<S>> {
        polar.dist =
            (self.dims.width() * self.dims.height()) /
            (
                pow(self.dims.width() * polar.t.sin(), 2) +
                pow(self.dims.height() * polar.t.cos(), 2)
            ).sqrt();
        polar
    }
}
