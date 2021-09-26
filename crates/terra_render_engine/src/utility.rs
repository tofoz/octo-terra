use super::pixel::*;

pub fn xy_to_index<T: std::ops::Mul<Output = T> + std::ops::Add<Output = T>>(
    x: T,
    y: T,
    width: T,
) -> T {
    x + y * width
}

pub fn index_to_xy<T: std::ops::Div<Output = T> + std::ops::Rem<Output = T> + Copy>(
    i: T,
    width: T,
) -> (T, T) {
    (i % width, i / width)
}

pub fn rgb_as_u32(col: &Pixel) -> u32 {
    ((col.a as u32) << 24) + ((col.r as u32) << 16) + ((col.g as u32) << 8) + (col.b as u32)
}

pub fn u32_as_argb(ucol: u32) -> Pixel {
    Pixel::new(
        (ucol << 24) as u8,
        (ucol << 16) as u8,
        (ucol << 8) as u8,
        ucol as u8,
    )
}
