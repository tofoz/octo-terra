extern crate num;
use super::*;
use std::{cmp::*, f32::NAN};

pub fn remap<
    T: std::ops::Sub<Output = T>
        + std::ops::Mul<Output = T>
        + std::ops::Div<Output = T>
        + std::ops::Add<Output = T>
        + Copy,
>(
    val: T,
    low1: T,
    high1: T,
    low2: T,
    high2: T,
) -> T {
    low2 + (val - low1) * (high2 - low2) / (high1 - low1)
}

pub fn remap_i32dp(val: i32, low1: i32, high1: i32, low2: i32, high2: i32) -> i32 {
    let hl2 = high2 - low2;
    let hl1 = high1 - low1;
    if hl1 == 0 || hl2 == 0 {
        low2
    } else {
        low2 + (val - low1) * (hl2) / (hl1)
    }
}

/// A value bounded by a minimum and a maximum
///
///  If input is less than min then this returns min.
///  If input is greater than max then this returns max.
///  Otherwise this returns input.
///
/// **Panics** in debug mode if `!(min <= max)`.
///
/// Original implementation from num-traits licensed as MIT
pub fn clamp<T: PartialOrd>(input: T, min: T, max: T) -> T {
    debug_assert!(min <= max, "min must be less than or equal to max");
    if input < min {
        min
    } else if input > max {
        max
    } else {
        input
    }
}

pub fn min3<T: Copy + PartialOrd>(a: T, b: T, c: T) -> T {
    let mut out_v = if a < b { a } else { b };
    out_v = if out_v < c { out_v } else { c };
    out_v
    //min(min(a, b), c)
}
pub fn max3<T: Copy + PartialOrd>(a: T, b: T, c: T) -> T {
    let mut out_v = if a > b { a } else { b };
    out_v = if out_v > c { out_v } else { c };
    out_v
    //max(a, max(b, c))
}

/// returns point of interest,
/// plane_n needs to be normalized.
pub fn vector_intersect_plane(
    plane_p: &Vec4,
    plane_n: &Vec4,
    line_start: &Vec4,
    line_end: &Vec4,
    t: &mut f32,
) -> Vec4 {
    let plane_n = plane_n.normalize();
    let plane_d = -plane_n.dot(*plane_p);
    let ad = line_start.dot(plane_n);
    let bd = line_end.dot(plane_n);
    *t = (-plane_d - ad) / (bd - ad);
    let line_start_to_end = *line_end - *line_start;
    let line_to_intersect = line_start_to_end * *t;
    *line_start + line_to_intersect
}
pub fn vector_intersect_plane_factor(
    plane_p: &Vec4,
    plane_n: &Vec4,
    line_start: &Vec4,
    line_end: &Vec4,
) -> f32 {
    let plane_n = plane_n.normalize();
    let plane_d = -plane_n.dot(*plane_p);
    let ad = line_start.dot(plane_n);
    let bd = line_end.dot(plane_n);
    (-plane_d - ad) / (bd - ad)
}

pub fn sign(x: f32) -> f32 {
    if x.is_sign_negative() {
        return -1.0;
    } else if x.is_sign_positive() {
        return 1.0;
    }
    if x == 0.0 {
        return 0.0;
    }

    NAN
}

pub fn sign_v3(v: Vec3) -> Vec3 {
    Vec3::new(sign(v.x), sign(v.y), sign(v.z))
}
