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

/// Gets the minimum value of all impute. Works on most types.
///
/// Each float impute must specify its type like so.
///``` rust
///
///     min!(2.5_f32, 6.1_f32, -3.71_f32);
///
///     fn foo(x: f64, y: f64, z: f64) {
///         min!(x, y, z);
///     }
///
///```
#[macro_export]
macro_rules! min {
    ($x:expr) => ($x);

    ($x:expr, $($y:expr), +) => (
    $x.min(min!($($y), +))
    )
}

/// Gets the maximum value of all impute. Works on most types.
///
/// Each float impute must specify its type like so.
///``` rust
///
///     max!(2.5_f32, 6.1_f32, -3.71_f32);
///
///     fn bar(x: f64, y: f64, z: f64) {
///         max!(x, y, z);
///     }
///
///```
#[macro_export]
macro_rules! max {
    ($x:expr) => ($x);

    ($x:expr, $($y:expr), +) => (
$x.max(max!($($y), +))
    )
}

/// prints the input expression fallowed by the result.
/// useful for debugging.
#[macro_export]
macro_rules! println_expression {
    ($e:expr) => {
        println!("{:?} = {:?}", stringify!($e), $e);
    };
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

mod testing {
    #[test]
    fn min_test() {
        assert_eq!(min!(2, -5, 1, 3, 0, -1), -5)
    }

    #[test]
    fn max_test() {
        assert_eq!(max!(2, -5, 1, 3, 0, -1), 3)
    }

    #[test]
    fn min_max_test() {
        assert_eq!(min!(15, 13, max!(-5, 2, 3, 7)), 7)
    }

    #[test]
    fn min_f32_test() {
        assert_eq!(
            min!(2.0_f32, -5.0_f32, 1.0_f32, 3.0_f32, 0.0_f32, -1.0_f32),
            -5.0_f32
        )
    }

    #[test]
    fn max_f32_test() {
        assert_eq!(
            max!(2.0_f32, -5.0_f32, 1.0_f32, 3.0_f32, 0.0_f32, -1.0_f32),
            3.0_f32
        )
    }
}
