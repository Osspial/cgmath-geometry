use crate::Dimensionality;
use std::iter::{ExactSizeIterator, FusedIterator, DoubleEndedIterator};
use super::GeoBox;

#[derive(Clone)]
pub struct NearestPointsIter<B: GeoBox> {
    pub(crate) b: B,
    pub(crate) point: <B::D as Dimensionality<B::Scalar>>::Point,
    /// Bitmask representing which sides of the box are closest. Each pair of two bits represent an
    /// axis, with the lowest bit in the pair representing the min side of the box and the highest
    /// bit representing the max side of the box.
    ///
    /// If this value == !0, return `Some(self.point)`.
    pub(crate) closest_sides: u16,
}

impl<B: GeoBox> NearestPointsIter<B> {
    fn get_with_side(&mut self, closest_side: u32) -> Option<<Self as Iterator>::Item> {
        if self.closest_sides == !0 {
            self.closest_sides = 0;
            Some(self.point)
        } else if self.closest_sides != 0 {
            self.closest_sides &= !(1 << closest_side);
            let closest_axis = (closest_side / 2) as usize;
            let mut point = self.point.clone();
            point[closest_axis] = match closest_side % 2 {
                0 => self.b.min()[closest_axis],
                1 => self.b.max()[closest_axis],
                _ => unreachable!()
            };
            Some(point)
        } else {
            None
        }
    }
}

impl<B: GeoBox> Iterator for NearestPointsIter<B> {
    type Item = <B::D as Dimensionality<B::Scalar>>::Point;
    fn next(&mut self) -> Option<Self::Item> {
        fn lowest_one_bit(n: u16) -> Option<u32> {
            for i in 0..16 {
                if (n >> i) & 1 == 1 {
                    return Some(i);
                }
            }
            None
        }

        let closest_side = lowest_one_bit(self.closest_sides).unwrap_or(0);
        self.get_with_side(closest_side)
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        if self.closest_sides == !0 {
            (1, Some(1))
        } else {
            let len = self.closest_sides.count_ones() as usize;
            (len, Some(len))
        }
    }
}

impl<B: GeoBox> DoubleEndedIterator for NearestPointsIter<B> {
    fn next_back(&mut self) -> Option<Self::Item> {
        fn highest_one_bit(n: u16) -> Option<u32> {
            for i in (0..16).rev() {
                if (n >> i) & 1 == 1 {
                    return Some(i);
                }
            }
            None
        }

        let closest_side = highest_one_bit(self.closest_sides).unwrap_or(0);
        self.get_with_side(closest_side)
    }
}
impl<B: GeoBox> ExactSizeIterator for NearestPointsIter<B> {}
impl<B: GeoBox> FusedIterator for NearestPointsIter<B> {}
