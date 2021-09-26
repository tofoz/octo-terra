extern crate prisma_maths;

use super::pixel::*;
use super::utility::*;
use prisma_maths::*;

//TODO: convert buffer from using pixels to using a generic type
// use u32 as the store value
// use u32 and convert u8 rgba back and forth, also can use f32
// add f32 mode for rgba chanel's or just a single r value ?
#[derive(Clone)]
pub struct RetroBitMap {
    pub buff_size: (usize, usize),
    pub buffer: Vec<Pixel>,
}

impl RetroBitMap {
    pub fn new(size: (usize, usize)) -> Self {
        RetroBitMap {
            buff_size: (size.0, size.1),
            buffer: vec![Pixel::new(0, 0, 0, 0); size.0 * size.1],
        }
    }

    pub fn clear_buffer(&mut self, color: &Pixel) {
        for c in self.buffer.iter_mut() {
            *c = *color;
        }
    }

    pub fn load_bitmap_from_image(img: &str) -> RetroBitMap {
        let dyn_img = image::open(img).unwrap();
        let img_size = image::image_dimensions(img).unwrap();
        let mut new_bit = RetroBitMap::new((img_size.0 as usize, img_size.1 as usize));

        let tt = dyn_img.as_rgba8().unwrap();

        for (i, p) in tt.pixels().enumerate() {
            //let (x, y) = index_to_xy(i, self.buff_size.0);

            new_bit.buffer[i] = Pixel::new(p[0], p[1], p[2], p[3]);
            //self.set_pixel(&Pixel::new(p[0], p[1], p[2], p[3]));
        }
        new_bit
    }

    // todo: add interpulason, etc
    // floor and cell and interpulat the diff
    // textur sampleing is off xy is only half
    pub fn sample_bit_map(&self, mut x: f32, mut y: f32) -> &Pixel {
        x %= 1.0;
        y %= 1.0;
        self.get_pixel(
            clamp(
                (x * 2.0) * self.buff_size.0 as f32,
                0.0,
                self.buff_size.0 as f32 - 1.,
            ) as u32,
            clamp(
                (y * 2.0) * self.buff_size.1 as f32,
                0.0,
                self.buff_size.1 as f32 - 1.,
            ) as u32,
        )
    }

    pub fn get_pixel(&self, x: u32, y: u32) -> &Pixel {
        let i = xy_to_index(x as usize, y as usize, self.buff_size.0);
        &self.buffer[i]
    }

    pub fn set_pixel(&mut self, x: i32, y: i32, col: &Pixel) {
        if let Some(mut p) = self
            .buffer
            .get_mut(xy_to_index(x, y, self.buff_size.0 as i32) as usize)
        {
            if x < self.buff_size.0 as i32 && x >= 0 {
                p.r = col.r;
                p.g = col.g;
                p.b = col.b;
                p.a = col.a;
            }
            //p = color;
        }
    }

    pub fn draw_rec(&mut self, rec: Rect<i32>, color: Pixel) {
        let p1 = rec.get_first_point();
        let p2 = rec.get_second_point();
        for i in p1.0..(p2.0 + 1) {
            for k in p1.1..(p2.1 + 1) {
                self.set_pixel(i, k, &color);
            }
        }
    }

    pub fn draw_bitmap(
        &mut self,
        img: &RetroBitMap,
        xi: i32,
        yi: i32,
        rec: (i32, i32, i32, i32),
        col0: Option<&Pixel>,
    ) {
        for x in rec.0..rec.0 + rec.2 {
            for y in rec.1..rec.1 + rec.3 {
                let yl = remap(y, rec.1, rec.1 + rec.3, rec.1 + rec.3, rec.1);
                let p = match col0 {
                    Some(c) => {
                        img.buffer[xy_to_index(x as usize, y as usize, img.buff_size.0)].blend(c)
                    }
                    None => img.buffer[xy_to_index(x as usize, y as usize, img.buff_size.0)],
                };

                if p.a > 15 {
                    self.set_pixel(xi + (x - rec.0), yi + (yl - rec.1), &p);
                }
            }
        }
    }

    // draws a bit map, can rotate, can scale, can squw rotate

    pub fn rot_copy() {}

    pub fn scale_copy(&self, scale: (f32, f32)) -> RetroBitMap {
        let mut nbm = RetroBitMap::new((
            (self.buff_size.0 as f32 * scale.0) as usize,
            (self.buff_size.1 as f32 * scale.1) as usize,
        ));
        for (i, p) in nbm.buffer.iter_mut().enumerate() {
            let (x, y) = index_to_xy(i, nbm.buff_size.0);
            let c = self.get_pixel((x as u32) / scale.0 as u32, (y as u32) / scale.1 as u32);
            *p = *c;
        }
        nbm
    }

    pub fn viwe_copy(&self, rec: (i32, i32, i32, i32)) -> RetroBitMap {
        let mut nbm = RetroBitMap::new((rec.2 as usize, rec.3 as usize));
        for x in rec.0..rec.0 + rec.2 {
            for y in rec.1..rec.1 + rec.3 {
                let yl = remap(y, rec.1, rec.1 + rec.3, rec.1 + rec.3, rec.1);
                let p = self
                    .buffer
                    .get(xy_to_index(x as usize, y as usize, self.buff_size.0));

                if let Some(p) = p {
                    if p.a > 15 {
                        nbm.set_pixel(x - rec.0, yl - rec.1, p);
                    }
                }
            }
        }
        nbm
    }

    pub fn draw_circle(&mut self, x: i32, y: i32, r: i32) {
        //clip bonding box and
        //use sdf to draw circle
    }

    pub fn draw_char(&mut self, img: &RetroBitMap, x: i32, y: i32, c: char, col: Option<&Pixel>) {
        let (ix, iy) = index_to_xy(c as i32, 16);
        self.draw_bitmap(img, x, y, (ix * 8, iy * 8, 8, 8), col);
    }

    pub fn draw_string(
        &mut self,
        img: &RetroBitMap,
        x: i32,
        y: i32,
        size: (u32, u32),
        string: &str,
        col: Option<&Pixel>,
    ) {
        for (i, s) in string.char_indices() {
            self.draw_char(img, x + (8 * i as i32), y, s, col)
        }
    }
}
