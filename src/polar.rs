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

use cgmath::*;

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Polar2<A: Angle> {
    /// Polar Angle (theta).
    pub t: A,
    pub dist: A::Unitless
}

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Polar3<A: Angle> {
    /// Polar Angle (theta), along XY plane.
    pub t: A,
    /// Azimuthal Angle (phi), from the Z axis.
    pub p: A,
    /// Distance from the origin.
    pub dist: A::Unitless
}

pub trait PolarSpace {
    type Scalar: BaseFloat;
    type Angle: Angle<Unitless=Self::Scalar>;
    type Array: AsRef<[Self::Angle]>;
    type Euclid: EuclideanSpace<Scalar=Self::Scalar>;

    fn angle_array(&self) -> &Self::Array;
    fn dist(&self) -> <Self::Angle as Angle>::Unitless;

    fn to_euclid(self) -> Self::Euclid;
    fn from_euclid(e: Self::Euclid) -> Self;
}

impl<A: Angle> PolarSpace for Polar2<A> {
    type Scalar= A::Unitless;
    type Angle = A;
    type Array = [Self::Angle; 1];
    type Euclid = Point2<A::Unitless>;

    #[inline(always)]
    fn angle_array(&self) -> &[Self::Angle; 1] {
        unsafe{ &*(&self.t as *const Self::Angle as *const [Self::Angle; 1] ) }
    }
    #[inline(always)]
    fn dist(&self) -> A::Unitless {
        self.dist
    }

    #[inline(always)]
    fn to_euclid(self) -> Point2<A::Unitless> {
        Point2 {
            x: self.dist * self.t.cos(),
            y: self.dist * self.t.sin()
        }
    }
    #[inline(always)]
    fn from_euclid(e: Point2<A::Unitless>) -> Polar2<A> {
        Polar2 {
            t: Self::Angle::atan(e.y/e.x),
            dist: e.to_vec().magnitude()
        }
    }
}

impl<A: Angle> PolarSpace for Polar3<A> {
    type Scalar = A::Unitless;
    type Angle = A;
    type Array = [A; 2];
    type Euclid = Point3<A::Unitless>;

    #[inline(always)]
    fn angle_array(&self) -> &[Self::Angle; 2] {
        unsafe{ &*(&self.t as *const Self::Angle as *const [Self::Angle; 2] ) }
    }
    #[inline(always)]
    fn dist(&self) -> A::Unitless {
        self.dist
    }

    #[inline(always)]
    fn to_euclid(self) -> Point3<A::Unitless> {
        Point3 {
            x: self.dist * self.t.cos() * self.p.sin(),
            y: self.dist * self.t.sin() * self.p.sin(),
            z: self.dist * self.p.cos()
        }
    }
    #[inline(always)]
    fn from_euclid(e: Point3<A::Unitless>) -> Polar3<A> {
        let dist = e.to_vec().magnitude();
        Polar3 {
            t: Self::Angle::atan(e.y/e.x),
            p: Self::Angle::acos(e.z/dist),
            dist
        }
    }
}
