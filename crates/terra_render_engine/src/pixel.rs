use std::ops::{Add, Div, Mul, Sub};

use prisma_maths::*;

#[derive(Copy, Clone)]
pub struct Pixel {
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub a: u8,
}

impl Add for Pixel {
    type Output = Pixel;

    fn add(self, rhs: Self) -> Self::Output {
        Pixel::new(
            self.r + rhs.r,
            self.g + rhs.g,
            self.b + rhs.b,
            self.a + rhs.a,
        )
    }
}
impl Sub for Pixel {
    type Output = Pixel;

    fn sub(self, rhs: Self) -> Self::Output {
        Pixel::new(
            self.r - rhs.r,
            self.g - rhs.g,
            self.b - rhs.b,
            self.a - rhs.a,
        )
    }
}
impl Mul for Pixel {
    type Output = Pixel;

    fn mul(self, rhs: Self) -> Self::Output {
        Pixel::new(
            self.r * rhs.r,
            self.g * rhs.g,
            self.b * rhs.b,
            self.a * rhs.a,
        )
    }
}
impl Div for Pixel {
    type Output = Pixel;

    fn div(self, rhs: Self) -> Self::Output {
        Pixel::new(
            self.r / rhs.r,
            self.g / rhs.g,
            self.b / rhs.b,
            self.a / rhs.a,
        )
    }
}

impl Pixel {
    pub fn new(r: u8, g: u8, b: u8, a: u8) -> Self {
        Pixel {
            r: r,
            g: g,
            b: b,
            a: a,
        }
    }
    pub fn get_col(&self) -> (u8, u8, u8, u8) {
        (self.r, self.g, self.b, self.a)
    }
    pub fn blend(&self, op: &Pixel) -> Pixel {
        Pixel::new(
            remap(op.a as i32, 0, 255, self.r as i32, op.r as i32) as u8,
            remap(op.a as i32, 0, 255, self.g as i32, op.g as i32) as u8,
            remap(op.a as i32, 0, 255, self.b as i32, op.b as i32) as u8,
            self.a,
        )
    }
    pub fn as_float(&self) -> (f32, f32, f32, f32) {
        (
            self.r as f32 / 255.,
            self.g as f32 / 255.,
            self.b as f32 / 255.,
            self.a as f32 / 255.,
        )
    }
}
