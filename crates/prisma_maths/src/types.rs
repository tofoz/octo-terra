extern crate num;

//Todo: create a sup folder/mod for types
//extern crate num_traits; // lets you do number ops on generics

//use num_traits::{pow, Signed};

use core::{
    cmp::*,
    f32,
    ops::{Add, Div, Mul, Sub},
};

use glam::Vec4;

struct foo {
    v: f32,
}

struct bar {
    x: f32,
    y: f32,
}

impl Add<bar> for foo {
    type Output = foo;

    fn add(self, rhs: bar) -> Self::Output {
        foo {
            v: self.v + rhs.x + rhs.y,
        }
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct Vector3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Vector3<T> {
    pub fn new(x: T, y: T, z: T) -> Self {
        Self { x, y, z }
    }
}

impl<T> Vector3<T>
where
    T: num::Signed + Copy + Mul<Output = T> + Add<Output = T>,
{
    pub fn dot(&self, other: &Vector3<T>) -> T {
        (self.x * other.x) + (self.y * other.y) + (self.z * other.z)
    }
    pub fn cross(&self, other: &Vector3<T>) -> Vector3<T> {
        Vector3::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }
    pub fn abs(&self) -> Vector3<T> {
        Vector3::new(num::abs(self.x), num::abs(self.y), num::abs(self.z))
    }
}

// functions that require floating point
impl<T> Vector3<T>
where
    T: num::Float + num::Signed,
{
    pub fn length(&self) -> T {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn normalize(&self) -> Vector3<T> {
        let len = self.length();
        if len > num::zero() {
            let inv_len = num::one::<T>() / len;
            return Vector3::new(self.x * inv_len, self.y * inv_len, self.z * inv_len);
        }
        Vector3::new(num::zero(), num::zero(), num::zero())
    }
}

impl<T> Add for Vector3<T>
where
    T: Add<Output = T>,
{
    type Output = Vector3<T>;

    fn add(self, rhs: Self) -> Self::Output {
        Vector3 {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}
impl<T> Sub for Vector3<T>
where
    T: Sub<Output = T>,
{
    type Output = Vector3<T>;

    fn sub(self, rhs: Self) -> Self::Output {
        Vector3 {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}
impl<T> Mul for Vector3<T>
where
    T: Mul<Output = T>,
{
    type Output = Vector3<T>;

    fn mul(self, rhs: Self) -> Self::Output {
        Vector3 {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
            z: self.z * rhs.z,
        }
    }
}
impl<T> Div for Vector3<T>
where
    T: Div<Output = T>,
{
    type Output = Vector3<T>;

    fn div(self, rhs: Self) -> Self::Output {
        Vector3 {
            x: self.x / rhs.x,
            y: self.y / rhs.y,
            z: self.z / rhs.z,
        }
    }
}

/// gets componet of vec4
pub fn v4_get(s: Vec4, index: i32) -> f32 {
    match index {
        0 => s.x,
        1 => s.y,
        2 => s.z,
        3 => s.w,
        _ => panic!("out of range!"),
    }
}
#[derive(Copy, Clone)]
pub struct Rect<T> {
    x: T,
    y: T,
    w: T,
    h: T,
}

impl<T: Mul<Output = T> + PartialOrd + Copy + Add<Output = T> + Ord + Sub<Output = T>> Rect<T> {
    pub fn new(x: T, y: T, w: T, h: T) -> Rect<T> {
        Rect { x, y, w, h }
    }

    pub fn get_wight(&self) -> T {
        self.w
    }
    pub fn get_hight(&self) -> T {
        self.h
    }
    pub fn set_wight(&mut self, w: T) {
        self.w = w;
    }
    pub fn set_hight(&mut self, h: T) {
        self.h = h;
    }
    pub fn get_first_point(&self) -> (T, T) {
        (self.x, self.y)
    }
    pub fn get_second_point(&self) -> (T, T) {
        (self.x + self.w, self.y + self.h)
    }
    pub fn set_first_point(&mut self, x: T, y: T) {
        self.x = x;
        self.y = y;
    }
    pub fn set_second_point(&mut self, x: T, y: T) {
        self.w = x - self.x;
        self.h = y - self.y;
    }
    ///gets a xy and b xy
    pub fn get_absolute(&self) -> (T, T, T, T) {
        (self.x, self.y, self.x + self.w, self.y + self.h)
    }
    pub fn set_absolute(&mut self, xa: T, ya: T, xb: T, yb: T) {
        self.x = xa;
        self.y = ya;
        self.w = xb - xa;
        self.h = yb - ya;
    }

    pub fn area(&self) -> T {
        self.w * self.h
    }

    pub fn can_hold(&self, other: &Rect<T>) -> bool {
        self.w > other.w && self.h > other.h
    }
    pub fn is_point_inside(&self, point: (T, T)) -> bool {
        point.0 > self.x
            && point.0 < self.x + self.w
            && point.1 > self.y
            && point.1 < self.y + self.h
    }
    pub fn is_overlap(&self, other: &Rect<T>) -> bool {
        !(self.x > other.x + other.w
        // R1 is right to R2
        || self.x + self.w < other.x
        // R1 is left to R2
        || self.y + self.h < other.y
        // R1 is above R2
        || self.y > other.y + other.h)
    }

    //
    pub fn get_overlap(&self, other: &Rect<T>) -> Option<Rect<T>> {
        let (x1, y1, x2, y2) = self.get_absolute();
        let (x3, y3, x4, y4) = other.get_absolute();
        let left = max(x1, x3); //x <-
        let right = min(x2, x4); // x ->
        let bot = max(y1, y3); // y \/
        let top = min(y2, y4); // y ^

        if left > right || bot > top {
            None
        } else {
            let mut rr = Rect::new(x1, x1, x1, x1);
            rr.set_absolute(left, bot, right, top);
            Some(rr)
        }
    }
}
