extern crate prisma_maths;

use crate::retro_bit_map;

use super::camera::*;
use super::mesh::*;
use super::pixel::*;
use super::ray_tracing::*;
use super::retro_bit_map::*;
use super::utility::*;
use prisma_maths::*;

use std::mem::swap;
use std::{
    cmp::*,
    mem,
    ops::{Add, Mul, Sub},
};

#[derive(Clone)]
pub struct RenderEngine {
    pub buff_size: (usize, usize),
    pub buffer: RetroBitMap, //TODO: replace with Retro Bit Map
    depth_buffer: Option<Vec<f32>>,
    vertices: Vec<Vertex>,
    aux_list: Vec<Vertex>,
}

impl RenderEngine {
    pub fn new(size: (usize, usize)) -> Self {
        RenderEngine {
            buff_size: (size.0, size.1),
            buffer: RetroBitMap::new(size),
            depth_buffer: None,
            vertices: Vec::new(),
            aux_list: Vec::new(),
        }
    }
    pub fn init_depth_buffer(&mut self) {
        self.depth_buffer = Some(vec![f32::INFINITY; self.buff_size.0 * self.buff_size.1])
    }

    pub fn clear_depth_buffer(&mut self) {
        match self.depth_buffer.as_mut() {
            Some(db) => {
                for f in db.iter_mut() {
                    *f = f32::INFINITY;
                }
            }
            None => {}
        }
    }

    pub fn get_depth_pixel(&self, x: u32, y: u32) -> f32 {
        let to = match self.depth_buffer.as_ref() {
            Some(db) => {
                let i = xy_to_index(x as usize, y as usize, self.buff_size.0);
                db.get(i)
            }
            None => None,
        };
        match to {
            Some(fv) => *fv,
            None => f32::INFINITY,
        }
    }

    pub fn set_depth_pixel(&mut self, x: i32, y: i32, dv: f32) {
        match self.depth_buffer.as_mut() {
            Some(db) => match db.get_mut(xy_to_index(x, y, self.buff_size.0 as i32) as usize) {
                Some(p) => {
                    if x < self.buff_size.0 as i32 && x >= 0 {
                        *p = dv;
                    }
                }
                None => {}
            },
            None => {}
        }
    }

    pub fn draw_rec(&mut self, rec: Rect<i32>, depth: f32, color: Pixel) {
        let p1 = rec.get_first_point();
        let p2 = rec.get_second_point();
        for i in p1.0..(p2.0 + 1) {
            for k in p1.1..(p2.1 + 1) {
                if self.get_depth_pixel(i as u32, k as u32) > depth {
                    self.buffer.set_pixel(i, k, &color);
                    self.set_depth_pixel(i, k, depth);
                }
            }
        }
    }

    pub fn draw_line(
        &mut self,
        xa: i32,
        ya: i32,
        za: f32,
        xb: i32,
        yb: i32,
        zb: f32,
        col0: Pixel,
        col1: Pixel,
    ) {
        let point_list = draw_line_pl(xa, ya, xb, yb);

        for (i, p) in point_list.0.iter().enumerate() {
            let s = if point_list.1 == 0 {
                remap(i, 0, point_list.0.len(), 0, 255)
            } else {
                remap(i as i32, 0, point_list.0.len() as i32, 255, 0) as usize
            };

            match self.depth_buffer.as_mut() {
                Some(db) => {
                    let ndv = remap(i as f32, 0.0, point_list.0.len() as f32, za, zb);
                    if self.get_depth_pixel(p.0 as u32, p.1 as u32) > ndv && ndv > 1.0 {
                        self.set_depth_pixel(p.0, p.1, ndv);

                        let nc = Pixel::new(col1.r, col1.g, col1.b, s as u8);

                        let c = col0.blend(&nc);
                        self.buffer.set_pixel(p.0, p.1, &c);
                    }
                }
                None => {
                    let nc = Pixel::new(col1.r, col1.g, col1.b, s as u8);

                    let c = col0.blend(&nc);
                    self.buffer.set_pixel(p.0, p.1, &c);
                }
            }
        }
    }

    pub fn dep_draw_triangle<T, V, F>(
        &mut self,
        v0: Vertex,
        v1: Vertex,
        v2: Vertex,
        points_inside: i32,
        shader_pipe: &mut SharderPipe<T, V, F>,
    ) where
        V: Fn(&mut Vertex, &mut T),
        F: Fn(Vertex, &mut T) -> Pixel,
    {
        let hw = self.buff_size.0 as f32 / 2.;
        let hh = self.buff_size.1 as f32 / 2.;
        if points_inside == 3 {
            let mut vv0 = v0;
            let mut vv1 = v1;
            let mut vv2 = v2;
            vv0.point.x = (vv0.point.x / vv0.point.w) * hw + hw;
            vv0.point.y = (vv0.point.y / vv0.point.w) * hh + hh;
            vv0.point.z /= vv0.point.w;

            vv1.point.x = (vv1.point.x / vv1.point.w) * hw + hw;
            vv1.point.y = (vv1.point.y / vv1.point.w) * hh + hh;
            vv1.point.z /= vv1.point.w;

            vv2.point.x = (vv2.point.x / vv2.point.w) * hw + hw;
            vv2.point.y = (vv2.point.y / vv2.point.w) * hh + hh;
            vv2.point.z /= vv2.point.w;

            self.draw_barycentric(&mut vv0, &mut vv1, &mut vv2, shader_pipe);
        } else {
            self.vertices.clear();
            self.aux_list.clear();
            self.vertices.push(v0);
            self.vertices.push(v1);
            self.vertices.push(v2);

            // only cliping the z axes for fps
            //# cliping is the hardest part so far!
            if
            /*RetroBitMap::clip_polygon_axis(&mut self.vertices, &mut self.aux_list, 0)
            && RetroBitMap::clip_polygon_axis(&mut self.vertices, &mut self.aux_list, 1)
            && */
            RenderEngine::clip_polygon_axis(&mut self.vertices, &mut self.aux_list, 2) {
                let mut initial_vertex = *self.vertices.get(0).unwrap();
                for i in 0..self.vertices.len() - 1 {
                    let mut vv0 = initial_vertex;
                    let mut vv1 = *self.vertices.get(i).unwrap();
                    let mut vv2 = *self.vertices.get(i + 1).unwrap();
                    // prospective divide
                    vv0.point.x = (vv0.point.x / vv0.point.w) * hw + hw;
                    vv0.point.y = (vv0.point.y / vv0.point.w) * hh + hh;
                    vv0.point.z /= vv0.point.w;

                    vv1.point.x = (vv1.point.x / vv1.point.w) * hw + hw;
                    vv1.point.y = (vv1.point.y / vv1.point.w) * hh + hh;
                    vv1.point.z /= vv1.point.w;

                    vv2.point.x = (vv2.point.x / vv2.point.w) * hw + hw;
                    vv2.point.y = (vv2.point.y / vv2.point.w) * hh + hh;
                    vv2.point.z /= vv2.point.w;

                    self.draw_barycentric(&mut vv0, &mut vv1, &mut vv2, shader_pipe);
                }
            }
        }
    }

    pub fn draw_triangle<T, V, F>(
        &mut self,
        mut v0: Vertex,
        mut v1: Vertex,
        mut v2: Vertex,
        shader_pipe: &mut SharderPipe<T, V, F>,
    ) where
        V: Fn(&mut Vertex, &mut T),
        F: Fn(Vertex, &mut T) -> Pixel,
    {
        let hw = self.buff_size.0 as f32 / 2.;
        let hh = self.buff_size.1 as f32 / 2.;

        let in_tri = Triangle::new((v0, v1, v2));

        let mut clipped: [Triangle; 2] = [in_tri; 2];

        let sp = clipped.split_at_mut(1); // may need to set to 1
        let mut n_clipped_triangles = RenderEngine::triangle_clip_against_plane(
            &Vec4::new(0., 0., 0.1, 0.),
            &Vec4::new(0., 0., 1.0, 0.).normalize(),
            &in_tri,
            &mut sp.0[0],
            &mut sp.1[0],
        );

        for n in 0..n_clipped_triangles as usize {
            v0 = clipped[n].p.0;
            v1 = clipped[n].p.1;
            v2 = clipped[n].p.2;

            // prospective divide
            v0.point.x = (v0.point.x / v0.point.w) * hw + hw;
            v0.point.y = (v0.point.y / v0.point.w) * hh + hh;
            v0.point.z /= v0.point.w;

            v1.point.x = (v1.point.x / v1.point.w) * hw + hw;
            v1.point.y = (v1.point.y / v1.point.w) * hh + hh;
            v1.point.z /= v1.point.w;

            v2.point.x = (v2.point.x / v2.point.w) * hw + hw;
            v2.point.y = (v2.point.y / v2.point.w) * hh + hh;
            v2.point.z /= v2.point.w;

            // check tri normal and flip if backwords
            let tri_norm = Vec3::cross(
                v1.point.xyz() - v0.point.xyz(),
                v2.point.xyz() - v0.point.xyz(),
            );

            let tri_facing = tri_norm.dot(vec3(0., 0., 1.0));
            if tri_facing <= 0.0 {
                mem::swap(&mut v0, &mut v1);
            }

            self.draw_barycentric(&mut v0, &mut v1, &mut v2, shader_pipe);
        }
    }

    //TODO: add materials / frag shadder, you can than add vertex lishting
    fn draw_barycentric<T, V, F>(
        &mut self,
        v0: &mut Vertex,
        v1: &mut Vertex,
        v2: &mut Vertex,
        shader_pipe: &mut SharderPipe<T, V, F>,
    ) where
        V: Fn(&mut Vertex, &mut T),
        F: Fn(Vertex, &mut T) -> Pixel,
    {
        //let mut v0 = *v0;
        //let mut v1 = *v1;
        //let mut v2 = *v2;

        // prospective divide here <-
        //todo: move out of here !!

        /*
        self.set_pixel(
            v0.point.x as i32 - 2,
            v0.point.y as i32,
            &Pixel::new(255, 00, 0, 255),
        );
        self.set_pixel(
            v1.point.x as i32 + 2,
            v1.point.y as i32,
            &Pixel::new(255, 255, 0, 255),
        );
        self.set_pixel(
            v2.point.x as i32,
            v2.point.y as i32 + 2,
            &Pixel::new(0, 00, 255, 255),
        );
        */

        //back face cull, wouldent draw regardles but try useing min and max insted to disregard try
        let tri_norm = Vec3::cross(
            v1.point.xyz() - v0.point.xyz(),
            v2.point.xyz() - v0.point.xyz(),
        );
        let tri_facing = tri_norm.dot(vec3(0., 0., 1.0));
        if tri_facing < 0.0 && shader_pipe.back_side_cull {
            return;
        } else if !shader_pipe.back_side_cull {
            if tri_facing < 0.0 {
                swap(v0, v2);
            }
        }

        // Compute triangle draw bounding box
        let mut min_x = min!(v0.point.x as i32, v1.point.x as i32, v2.point.x as i32);
        let mut min_y = min!(v0.point.y as i32, v1.point.y as i32, v2.point.y as i32);
        let mut max_x = max!(v0.point.x as i32, v1.point.x as i32, v2.point.x as i32);
        let mut max_y = max!(v0.point.y as i32, v1.point.y as i32, v2.point.y as i32);

        // Draw Clip against screen bounds

        min_x = max!(min_x, 0);
        min_y = max!(min_y, 0);
        max_x = min!(max_x, self.buff_size.0 as i32 - 1);
        max_y = min!(max_y, self.buff_size.1 as i32 - 1);

        let e0 = EdgeEquation::new(&v0.point.xyz(), &v1.point.xyz());
        let e1 = EdgeEquation::new(&v1.point.xyz(), &v2.point.xyz());
        let e2 = EdgeEquation::new(&v2.point.xyz(), &v0.point.xyz());
        let area = edge_function(
            &(v0.point.x, v0.point.y),
            &(v1.point.x, v1.point.y),
            &(v2.point.x, v2.point.y),
        );

        let r = ParameterEquation::new(v2.color.x, v0.color.x, v1.color.x, &e0, &e1, &e2, area);
        let g = ParameterEquation::new(v2.color.y, v0.color.y, v1.color.y, &e0, &e1, &e2, area);
        let b = ParameterEquation::new(v2.color.z, v0.color.z, v1.color.z, &e0, &e1, &e2, area);
        let uvx = ParameterEquation::new(v2.uv.x, v0.uv.x, v1.uv.x, &e0, &e1, &e2, area);
        let uvy = ParameterEquation::new(v2.uv.y, v0.uv.y, v1.uv.y, &e0, &e1, &e2, area);
        //let nz = ParameterEquation::new(v2.normal.z, v0.normal.z, v1.normal.z, &e0, &e1, &e2, area);
        let nx = ParameterEquation::new(v2.normal.x, v0.normal.x, v1.normal.x, &e0, &e1, &e2, area);
        let ny = ParameterEquation::new(v2.normal.y, v0.normal.y, v1.normal.y, &e0, &e1, &e2, area);
        let nz = ParameterEquation::new(v2.normal.z, v0.normal.z, v1.normal.z, &e0, &e1, &e2, area);

        let dp = ParameterEquation::new(v2.point.z, v0.point.z, v1.point.z, &e0, &e1, &e2, area);

        // todo: make multi tharaded with rayden?
        for y in min_y..max_y + 1 {
            for x in min_x..max_x + 1 {
                let p = (x as f32 + 0.5, y as f32 + 0.5);
                // get v0 val from the of=peset side langth?

                if e0.test_edge(p.0, p.1) && e1.test_edge(p.0, p.1) && e2.test_edge(p.0, p.1) {
                    /// pixel shafer ->
                    match self.depth_buffer.as_ref() {
                        Some(_) => {
                            let cdp = dp.evaluate(p.0, p.1);

                            if self.get_depth_pixel(x as u32, y as u32) > cdp
                            /*&& cdp < 0.5 */
                            {
                                let v_norm = Vec4::new(
                                    nx.evaluate(p.0, p.1),
                                    ny.evaluate(p.0, p.1),
                                    nz.evaluate(p.0, p.1),
                                    1.0,
                                );
                                let v_color = Vec4::new(
                                    r.evaluate(p.0, p.1),
                                    g.evaluate(p.0, p.1),
                                    b.evaluate(p.0, p.1),
                                    1.0,
                                );
                                let v_uv =
                                    Vec2::new(uvx.evaluate(p.0, p.1), uvy.evaluate(p.0, p.1));

                                // todo add uv as peramater
                                let vert = Vertex {
                                    point: Vec4::zero(),
                                    normal: v_norm.xyz(),
                                    color: v_color.xyz(),
                                    uv: v_uv,
                                };
                                let col = shader_pipe.run_frag_shader(vert);
                                if col.a >= (0.5 * 255.0) as u8 {
                                    self.buffer.set_pixel(x, y, &col);
                                    self.set_depth_pixel(x, y, cdp);
                                }
                            }
                        }
                        None => {
                            let norm = Vec4::new(
                                nx.evaluate(p.0, p.1),
                                ny.evaluate(p.0, p.1),
                                nz.evaluate(p.0, p.1),
                                1.0,
                            );
                            let light_level = clamp(
                                norm.xyz().dot(vec3(1.0, 1.0, 1.0).normalize()) * 0.6 + 0.4,
                                0.0,
                                1.0,
                            ) * 2.0;

                            self.buffer.set_pixel(
                                x,
                                y,
                                &Pixel::new(
                                    (r.evaluate(p.0, p.1) * light_level * 255.0) as u8,
                                    (g.evaluate(p.0, p.1) * light_level * 255.0) as u8,
                                    (b.evaluate(p.0, p.1) * light_level * 255.0) as u8,
                                    255,
                                ),
                            );
                        }
                    }
                }
            }
        }
    }

    // todo: create a shader_pipe insted of passing in vertex and frag shaders derectly

    pub fn draw_mesh<V, F, T>(
        &mut self,
        mesh: &Mesh,
        shader_pipe: &mut SharderPipe<T, V, F>,
        camera: &mut Camera,
        mode: &MeshRenderMode,
    ) where
        V: Fn(&mut Vertex, &mut T),
        F: Fn(Vertex, &mut T) -> Pixel,
    {
        //let nx = remap(0., 0., self.buff_size.0 as f32, -1., 1.);
        //let ny = remap(0., 0., self.buff_size.1 as f32, -1., 1.);

        // half size of screen is used for normalizason
        // normalized point [-1, 1] * half size + half size
        //TODO: move hw and hh out as varbs of render
        //todo: make cam object and add a internal on in render
        //

        let hw = self.buff_size.0 as f32 / 2.;
        let hh = self.buff_size.1 as f32 / 2.;

        let pm = camera.get_perspective_mat4() * camera.get_camera_to_world();

        //# object to world space -> world space to cam space -> cam space to nds/screen space

        match mode {
            MeshRenderMode::VoxelPoint => { /*
                 for i in 0..(m.v_index.len()) {
                     let mut a = *m.vertices.get(m.v_index[i]).unwrap();

                     a.point = m.transform * a.point;

                     //let pm = Mat4::perspective_rh_gl(fov.to_radians(), ar, 1.0, 100.0);
                     //let new_pm = pm * *cam_transform;

                     // prospective projection
                     a.point = pm * a.point;

                     if a.point.z < -a.point.w || a.point.z > a.point.w {
                         continue;
                     }
                     if a.point.y < -a.point.w || a.point.y > a.point.w {
                         continue;
                     }
                     if a.point.x < -a.point.w || a.point.x > a.point.w {
                         continue;
                     }

                     a.point.x = (a.point.x / a.point.w) * hw + hw;
                     a.point.y = (a.point.y / a.point.w) * hh + hh;

                     let mut rec = Rect::new(
                         a.point.x as i32,
                         a.point.y as i32,
                         (50.0 / a.point.w) as i32,
                         (50.0 / a.point.w) as i32,
                     );
                     let tt = rec.get_first_point();
                     let mp = (tt.0 - (rec.get_wight() / 2), tt.1 - (rec.get_hight() / 2));
                     rec.set_first_point(mp.0, mp.1);
                     let cdp = a.point.z / a.point.w;

                     let light_level =
                         clamp(a.normal.dot(vec3(0.0, 0.0, 1.0)) * 0.5 + 0.5, 0.0, 1.0) * 2.0;
                     self.draw_rec(
                         rec,
                         cdp,
                         Pixel::new(
                             (a.color.x * light_level * 255.0) as u8,
                             (a.color.y * light_level * 255.0) as u8,
                             (a.color.z * light_level * 255.0) as u8,
                             200,
                         ),
                     )
                 }
                 */
            }
            MeshRenderMode::Dots => {
                /*
                for i in 0..(m.v_index.len()) {
                    let mut a = *m.vertices.get(m.v_index[i]).unwrap();
                    let mut a = vec4(a.x, a.y, a.z, 1.0);
                    a = m.transform * a;

                    //let pm = Mat4::perspective_rh_gl(fov.to_radians(), ar, 1.0, 100.0);
                    //let new_pm = pm * *cam_transform;

                    // prospective projection
                    a = pm * a;
                    a.x = (a.x / a.w) * hw + hw;
                    a.y = (a.y / a.w) * hh + hh;

                    let cdp = a.z;
                    let zcol = 1. / (cdp * 0.075);
                    if self.get_depth_pixel(a.x as u32, a.y as u32) > cdp && cdp > 0.4 {
                        self.set_pixel(
                            (a.x) as i32,
                            (a.y) as i32,
                            &Pixel::new(
                                (m.color[m.v_index[i]].x * zcol * 255.0) as u8,
                                (m.color[m.v_index[i]].y * zcol * 255.0) as u8,
                                (m.color[m.v_index[i]].z * zcol * 255.0) as u8,
                                200,
                            ),
                        );
                        self.set_depth_pixel(a.x as i32, a.y as i32, cdp);
                    }
                }
                */
            }
            MeshRenderMode::Lines => {
                /*
                for i in 0..(m.v_index.len() / 2) {
                    let mut a = *m.vertices.get(m.v_index[(i * 2)]).unwrap();
                    let mut a = vec4(a.x, a.y, a.z, 1.0);
                    a = m.transform * a;

                    let mut b = *m.vertices.get(m.v_index[(i * 2) + 1]).unwrap();
                    let mut b = vec4(b.x, b.y, b.z, 1.0);
                    b = m.transform * b;

                    a.x = (a.x / (ar * a.z * tanhfov)) * hw + hw;
                    a.y = (a.y / (a.z * tanhfov)) * hh + hh;

                    b.x = (b.x / (ar * b.z * tanhfov)) * hw + hw;
                    b.y = (b.y / (b.z * tanhfov)) * hh + hh;

                    let ap = a.xyz();
                    let bp = b.xyz();

                    self.draw_line(
                        a.x as i32,
                        a.y as i32,
                        a.z,
                        b.x as i32,
                        b.y as i32,
                        b.z,
                        Pixel::new(255, 10, 50, 200),
                        Pixel::new(50, 255, 0, 200),
                    );
                }
                */
            }
            MeshRenderMode::WireFrame => {
                /*
                for i in 0..(m.v_index.len() / 3) {
                    let mut a = *m.vertices.get(m.v_index[(i * 3)]).unwrap();
                    let mut a = vec4(a.x, a.y, a.z, 1.0);
                    a = m.transform * a;

                    let mut b = *m.vertices.get(m.v_index[(i * 3) + 1]).unwrap();
                    let mut b = vec4(b.x, b.y, b.z, 1.0);
                    b = m.transform * b;

                    let mut c = *m.vertices.get(m.v_index[(i * 3) + 2]).unwrap();
                    let mut c = vec4(c.x, c.y, c.z, 1.0);
                    c = m.transform * c;

                    let ap = a.xyz();
                    let bp = b.xyz();
                    let cp = c.xyz();

                    a = pm * a;
                    b = pm * b;
                    c = pm * c;

                    // clip space here <--

                    // prospective divide p.xyz / p.w and remape to screen space
                    a.x = (a.x / a.w) * hw + hw;
                    a.y = (a.y / a.w) * hh + hh;

                    b.x = (b.x / b.w) * hw + hw;
                    b.y = (b.y / b.w) * hh + hh;

                    c.x = (c.x / c.w) * hw + hw;
                    c.y = (c.y / c.w) * hh + hh;

                    a.z = a.z / a.w;
                    b.z = b.z / b.w;
                    c.z = c.z / c.w;

                    let min_tri_z = min3(a.z, b.z, c.z);
                    let max_tri_z = max3(a.z, b.z, c.z);

                    //TODO: timp hack, need to add tri clipping

                    if min_tri_z <= -1.0 {
                        continue;
                    }

                    if max_tri_z >= 1.0 {
                        continue;
                    }

                    let tri_norm = Vec3::cross(bp - ap, cp - ap);

                    let tri_facing = tri_norm.dot(vec3(0., 0., 1.));

                    if tri_facing <= 0. {
                        self.draw_line(
                            a.x as i32,
                            a.y as i32,
                            a.z,
                            b.x as i32,
                            b.y as i32,
                            b.z,
                            Pixel::new(255, 10, 50, 200),
                            Pixel::new(50, 255, 0, 200),
                        );

                        self.draw_line(
                            b.x as i32,
                            b.y as i32,
                            b.x,
                            c.x as i32,
                            c.y as i32,
                            c.z,
                            Pixel::new(0, 255, 50, 200),
                            Pixel::new(50, 10, 255, 200),
                        );

                        self.draw_line(
                            c.x as i32,
                            c.y as i32,
                            c.z,
                            a.x as i32,
                            a.y as i32,
                            a.z,
                            Pixel::new(255, 10, 50, 200),
                            Pixel::new(50, 10, 255, 200),
                        );
                    }
                }
                */
            }
            MeshRenderMode::Tri => {
                // todo: make multi threaded
                if true {
                    for i in 0..(mesh.v_index.len() / 3) {
                        //todo: 3d vertex gets converted to 4d point here
                        let mut a = *mesh.vertices.get(mesh.v_index[(i * 3)]).unwrap();
                        let mut b = *mesh.vertices.get(mesh.v_index[(i * 3) + 1]).unwrap();
                        let mut c = *mesh.vertices.get(mesh.v_index[(i * 3) + 2]).unwrap();

                        // transform model points from local/object-space to world-space
                        // todo add vertex shader here
                        // its the vertex shaders resposabilaty to do porjecon mult?

                        //# recomputing a vertex of shard tri points
                        //# have a buffer for of a meshs vertices, v_buff, that we para loop over befor we clip?
                        //# buffer will have a known size of the mesh vert size
                        //? is it worth doing this or just parlal looping over all triangles?

                        shader_pipe.run_vertex_shader(&mut a);
                        shader_pipe.run_vertex_shader(&mut b);
                        shader_pipe.run_vertex_shader(&mut c);

                        b.point = mesh.transform * b.point;
                        a.point = mesh.transform * a.point;
                        c.point = mesh.transform * c.point;

                        b.normal = (mesh.transform
                            * Vec4::new(b.normal.x, b.normal.y, b.normal.z, 0.))
                        .xyz();
                        a.normal = (mesh.transform
                            * Vec4::new(a.normal.x, a.normal.y, a.normal.z, 0.))
                        .xyz();
                        c.normal = (mesh.transform
                            * Vec4::new(c.normal.x, c.normal.y, c.normal.z, 0.))
                        .xyz();

                        // V shader end <--------------------------------------------->

                        // prospective mult, now in 3d homogeneous clipping space
                        a.point = pm * a.point;
                        b.point = pm * b.point;
                        c.point = pm * c.point;

                        // if all points are outside of volume skip tri
                        if !((-a.point.w < a.point.z && a.point.z < a.point.w)
                            || (-b.point.w < b.point.z && b.point.z < b.point.w)
                            || (-c.point.w < c.point.z && c.point.z < c.point.w))
                        {
                            continue;
                        }

                        if !((-a.point.w < a.point.y && a.point.y < a.point.w)
                            || (-b.point.w < b.point.y && b.point.y < b.point.w)
                            || (-c.point.w < c.point.y && c.point.y < c.point.w))
                        {
                            continue;
                        }

                        if !((-a.point.w < a.point.x && a.point.x < a.point.w)
                            || (-b.point.w < b.point.x && b.point.x < b.point.w)
                            || (-c.point.w < c.point.x && c.point.x < c.point.w))
                        {
                            continue;
                        }

                        // count points so we only run what we need to throw the clipper
                        // no real fps improvement for chekiang if 0 points are inside the volume and bug were triangle is culled some times
                        let mut points_inside = 0;
                        if (-a.point.w < a.point.z && a.point.z < a.point.w)
                            && (-a.point.w < a.point.y && a.point.y < a.point.w)
                            && (-a.point.w < a.point.x && a.point.x < a.point.w)
                        {
                            points_inside += 1;
                        }
                        if (-b.point.w < b.point.z && b.point.z < b.point.w)
                            && (-b.point.w < b.point.y && b.point.y < b.point.w)
                            && (-b.point.w < b.point.x && b.point.x < b.point.w)
                        {
                            points_inside += 1;
                        }
                        if (-c.point.w < c.point.z && c.point.z < c.point.w)
                            && (-c.point.w < c.point.y && c.point.y < c.point.w)
                            && (-c.point.w < c.point.x && c.point.x < c.point.w)
                        {
                            points_inside += 1;
                        }

                        match 2 {
                            // cull triangle if one of its point is out of z bounds
                            0 => {
                                if (-a.point.w > a.point.z || a.point.z > a.point.w)
                                    || (-b.point.w > b.point.z || b.point.z > b.point.w)
                                    || (-c.point.w > c.point.z || c.point.z > c.point.w)
                                {
                                    continue;
                                }

                                // prospective divide
                                a.point.x = (a.point.x / a.point.w) * hw + hw;
                                a.point.y = (a.point.y / a.point.w) * hh + hh;

                                b.point.x = (b.point.x / b.point.w) * hw + hw;
                                b.point.y = (b.point.y / b.point.w) * hh + hh;

                                c.point.x = (c.point.x / c.point.w) * hw + hw;
                                c.point.y = (c.point.y / c.point.w) * hh + hh;

                                //self.draw_barycentric(&mut a, &mut b, &mut c);
                            }
                            // new triangle clipping
                            1 => {
                                self.draw_triangle(a, b, c, shader_pipe);
                            }
                            // old triangle clipping
                            2 => {
                                self.dep_draw_triangle(a, b, c, points_inside, shader_pipe);
                            }

                            _ => {}
                        }
                    }
                } else {
                }
            }
        }
    }

    /// plane_n needs to be normalized
    fn triangle_clip_against_plane(
        plane_p: &Vec4,
        plane_n: &Vec4,
        in_tri: &Triangle,
        out_tri1: &mut Triangle,
        out_tri2: &mut Triangle,
    ) -> i32 {
        let dist = |p: &Vec4| {
            return plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - plane_n.dot(*plane_p);
        };

        let mut inside_points = [in_tri.p.0; 3];
        let mut n_inside_point_count = 0;
        let mut outside_points = [in_tri.p.0; 3];
        let mut n_out_side_point_count = 0;

        let d0 = dist(*&&in_tri.p.0.point);
        let d1 = dist(*&&in_tri.p.1.point);
        let d2 = dist(*&&in_tri.p.2.point);

        if d0 >= 0.0 {
            inside_points[n_inside_point_count] = in_tri.p.0;
            n_inside_point_count += 1;
        } else {
            outside_points[n_out_side_point_count] = in_tri.p.0;
            n_out_side_point_count += 1;
        }
        if d1 >= 0.0 {
            inside_points[n_inside_point_count] = in_tri.p.1;
            n_inside_point_count += 1;
        } else {
            outside_points[n_out_side_point_count] = in_tri.p.1;
            n_out_side_point_count += 1;
        }
        if d2 >= 0.0 {
            inside_points[n_inside_point_count] = in_tri.p.2;
            n_inside_point_count += 1;
        } else {
            outside_points[n_out_side_point_count] = in_tri.p.2;
            n_out_side_point_count += 1;
        }

        if n_inside_point_count == 0 {
            return 0;
        }

        if n_inside_point_count == 3 {
            *out_tri1 = *in_tri;
            return 1;
        }

        // here!!! tri out 1 point order is importent
        if n_inside_point_count == 1 && n_out_side_point_count == 2 {
            out_tri1.p.0 = inside_points[0];

            let mut t = vector_intersect_plane_factor(
                &plane_p,
                &plane_n,
                &inside_points[0].point,
                &outside_points[0].point,
            );
            out_tri1.p.1 = inside_points[0].leap(outside_points[0], t);
            t = vector_intersect_plane_factor(
                &plane_p,
                &plane_n,
                &inside_points[0].point,
                &outside_points[1].point,
            );
            out_tri1.p.2 = inside_points[0].leap(outside_points[1], t);

            return 1;
        }

        if n_inside_point_count == 2 && n_out_side_point_count == 1 {
            out_tri1.p.0 = inside_points[0];
            out_tri1.p.1 = inside_points[1];

            let mut t = vector_intersect_plane_factor(
                &plane_p,
                &plane_n,
                &inside_points[0].point,
                &outside_points[0].point,
            );
            out_tri1.p.2 = inside_points[0].leap(outside_points[0], t);

            // if glish its here
            out_tri2.p.0 = inside_points[1];
            out_tri2.p.1 = out_tri1.p.2;
            t = vector_intersect_plane_factor(
                &plane_p,
                &plane_n,
                &inside_points[1].point,
                &outside_points[0].point,
            );
            out_tri2.p.2 = inside_points[1].leap(outside_points[0], t);

            return 2;
        }

        return -1;
    }

    // old clip
    fn clip_polygon_axis(
        vertices: &mut Vec<Vertex>,
        aux_list: &mut Vec<Vertex>,
        component_index: i32,
    ) -> bool {
        RenderEngine::clip_polygon_component(vertices, component_index, 1.0, aux_list);
        vertices.clear();

        if aux_list.is_empty() {
            return false;
        }

        RenderEngine::clip_polygon_component(aux_list, component_index, -1.0, vertices);
        aux_list.clear();

        return !vertices.is_empty();
    }

    fn clip_polygon_component(
        vertices: &Vec<Vertex>,
        component_index: i32,
        component_factor: f32,
        result: &mut Vec<Vertex>,
    ) {
        let mut prev_vertex = vertices.get(vertices.len() - 1).unwrap();
        let mut prev_component = v4_get(prev_vertex.point, component_index) * component_factor;
        let mut prev_inside = prev_component <= prev_vertex.point.w;

        for current_vertex in vertices.iter() {
            let current_component =
                v4_get(current_vertex.point, component_index) * component_factor;
            let current_inside = current_component <= current_vertex.point.w;

            if current_inside ^ prev_inside {
                let lerp_amt = (prev_vertex.point.w - prev_component)
                    / ((prev_vertex.point.w - prev_component)
                        - (current_vertex.point.w - current_component));
                result.push(prev_vertex.leap(*current_vertex, lerp_amt));
            }

            if current_inside {
                result.push(*current_vertex);
            }

            prev_vertex = current_vertex;
            prev_component = current_component;
            prev_inside = current_inside;
        }
    }

    pub fn voxel_ray_trace(
        &mut self,
        voxel_texture: &VoxelTexture,
        model_mat: &Mat4,
        camera: &mut Camera,
    ) {
        let steps: i32 = 256;
        let voxel_per_meter: f32 = 8.;
        let epsy: f32 = 0.0001;

        //if you have more then one volume you want to trace you would want to check if the volume you tracing is overlapping another
        // A: if overlapping trace both or more and use a z buff to determine visibility
        // B: if overlapping trace both or more at the same time?

        // ray per pixel loop
        for i in 0..self.buff_size.0 * self.buff_size.1 {
            let (x, y) = self::index_to_xy(i, self.buff_size.0);
            // gives a value(uv) normalized between [-1, 1]
            let uv = (Vec2::new(x as f32, y as f32) + Vec2::new(0.5, 0.5))
                / Vec2::new(self.buff_size.0 as f32, self.buff_size.1 as f32)
                * 2.0
                - Vec2::splat(1.0);

            let ray = camera.get_camera_ray(uv);

            let row: Vec3 = ray.origin;
            let rdw: Vec3 = ray.direction;
            let txx: Mat4 = *model_mat;
            let txi: Mat4 = model_mat.inverse();
            let rad: Vec3 = Vec3::new(
                voxel_texture.size.0 as f32 / voxel_per_meter / 2.0 - 0.0001,
                voxel_texture.size.1 as f32 / voxel_per_meter / 2.0 - 0.0001,
                voxel_texture.size.2 as f32 / voxel_per_meter / 2.0 - 0.0001,
            );
            let mut ot = Vec2::new(0., 0.);
            let mut on = Vec3::new(0., 0., 0.);
            let mut ou = Vec2::new(0., 0.);
            let mut of = 0;

            if box_intersect(
                &row, &rdw, &txx, &txi, &rad, &mut ot, &mut on, &mut ou, &mut of,
            ) {
                //self.buffer
                //    .set_pixel(x as i32, y as i32, &Pixel::new(255, 10, 155, 255));
                //continue;
                let mut rd = ray.direction;
                let mut ro = ray.origin + (ray.direction * ot.x.max(0.0));

                if self.get_depth_pixel(x as u32, y as u32) >= ot.x.max(0.0) {
                    //    break;
                }

                rd = (*model_mat * rd.extend(0.0)).xyz();

                ro = (*model_mat * ro.extend(1.0)).xyz() * voxel_per_meter;
                ro = ro
                    + (Vec3::new(
                        voxel_texture.size.0 as f32,
                        voxel_texture.size.1 as f32,
                        voxel_texture.size.2 as f32,
                    ) / 2.0);

                let mut nray = Ray::create_ray(ro, rd);

                // texture look up
                // 0 is just ray marching
                // 1 is branchless voxel traversal, weird axis related bug
                // 2 Branchless DDA Voxel Raycasting
                // https://www.shadertoy.com/view/4dX3zl
                //TODO: look into multi rez mip map levels to traverse to save performance
                match 2 {
                    0 => {
                        let mut t_step = 1.0 as f32;
                        for i in 0..steps {
                            let p = nray.origin + nray.direction * t_step; //;
                            t_step += 1.0;
                            // flooring so it samples correct unit, fixed out of bounds sample bug
                            let vv = match voxel_texture.get_voxel_index(
                                p.x.floor() as i32,
                                p.y.floor() as i32,
                                p.z.floor() as i32,
                            ) {
                                Some(voxi) => voxi,
                                None => {
                                    break;
                                }
                            };
                            let vl = voxel_texture.map.get(vv as usize); //# <--

                            match vl {
                                Some(v_lup) => {
                                    if *v_lup != 0 {
                                        let col =
                                            *voxel_texture.palette.get(*v_lup as usize).unwrap()
                                                as i32;
                                        let rr = (col) as u8;
                                        let gg = (col >> 8) as u8;
                                        let bb = (col >> 16) as u8;

                                        self.buffer.set_pixel(
                                            x as i32,
                                            y as i32,
                                            &Pixel::new(rr, gg, bb, 255),
                                        );
                                        break;
                                    } else {
                                        continue;
                                    }
                                }

                                None => {
                                    break;
                                }
                            }
                        }
                    }
                    1 => {
                        let mut origin = nray.origin;
                        let direction = nray.direction;

                        // bug somewere here <-
                        let mut point = origin.floor();
                        let step =
                            Vec3::new(sign(direction.x), sign(direction.y), sign(direction.z));
                        let t_delta = 1.0 / direction.abs().max(Vec3::splat(epsy));

                        // origanly splat(1.0) - ,but 0.5 seams to make the bug less visabul.
                        let mut t_max: Vec3 = (Vec3::splat(1.0)
                            - Vec3::new(
                                (origin.x * step.x).fract(),
                                (origin.y * step.y).fract(),
                                (origin.z * step.z).fract(),
                            ))
                            * t_delta;

                        for _ in 0..steps {
                            let voxid = voxel_texture.get_voxel_index(
                                point.x.floor() as i32,
                                point.y.floor() as i32,
                                point.z.floor() as i32,
                            );

                            match voxid {
                                Some(vi) => {
                                    let vox = *voxel_texture.map.get(vi as usize).unwrap();
                                    if vox != 0 {
                                        //return Some((point.x as f32, point.y as f32, point.z as f32));

                                        let col = *voxel_texture.palette.get(vox as usize).unwrap()
                                            as i32;
                                        let rr = (col) as u8;
                                        let gg = (col >> 8) as u8;
                                        let bb = (col >> 16) as u8;

                                        self.buffer.set_pixel(
                                            x as i32,
                                            y as i32,
                                            &Pixel::new(rr, gg, bb, 255),
                                        );
                                        break;
                                    } else {
                                    }
                                }
                                None => {
                                    break;
                                }
                            }

                            // performance hit here <-

                            let ab = (t_max.x < t_max.y, t_max.y < t_max.z, t_max.z < t_max.x);
                            let bb = (t_max.x <= t_max.z, t_max.y <= t_max.x, t_max.z <= t_max.y);

                            let s = (ab.0 && bb.0, ab.1 && bb.1, ab.2 && bb.2);

                            if s.0 {
                                point.x += step.x;
                                t_max.x += t_delta.x;
                            }
                            if s.1 {
                                point.y += step.y;
                                t_max.y += t_delta.y;
                            }
                            if s.2 {
                                point.z += step.z;
                                t_max.z += t_delta.z;
                            }
                        }
                    }
                    2 => {
                        let mut map_pos = nray.origin.floor();

                        let delta_dist = (nray.direction.length() / nray.direction).abs();

                        let ray_step = sign_v3(nray.direction);

                        let mut side_dist = (sign_v3(nray.direction) * (map_pos - nray.origin)
                            + ((sign_v3(nray.direction) * 0.5) + Vec3::splat(0.5)))
                            * delta_dist;

                        let mut mask = BVec3::new(false, false, false);
                        for _ in 0..steps {
                            let voxid = voxel_texture.get_voxel_index(
                                map_pos.x as i32,
                                map_pos.y as i32,
                                map_pos.z as i32,
                            );

                            match voxid {
                                Some(vi) => {
                                    let vox = voxel_texture.map.get(vi as usize); // map

                                    match vox {
                                        Some(vox) => {
                                            if *vox != 0 {
                                                // vox != 0
                                                //return Some((point.x as f32, point.y as f32, point.z as f32));

                                                let col = *voxel_texture
                                                    .palette
                                                    .get(*vox as usize)
                                                    .unwrap()
                                                    as i32;
                                                let mut rr = (col) as u8;
                                                let mut gg = (col >> 8) as u8;
                                                let mut bb = (col >> 16) as u8;

                                                if true {
                                                    let n = Vec3::select(
                                                        mask,
                                                        Vec3::new(1.0, 1.0, 1.0),
                                                        Vec3::new(-1.0, -1.0, -1.0),
                                                    )
                                                    .normalize();
                                                    let ll = (n
                                                        .dot(Vec3::new(1.0, 1.0, 0.2).normalize())
                                                        * 0.5
                                                        + 0.5)
                                                        * 22.;
                                                    rr = rr.saturating_add(ll as u8);
                                                    gg = gg.saturating_add(ll as u8);
                                                    bb = bb.saturating_add(ll as u8);
                                                }

                                                self.buffer.set_pixel(
                                                    x as i32,
                                                    y as i32,
                                                    &Pixel::new(rr, gg, bb, 255),
                                                );

                                                // todo: to get the depth we need the distance of the ray from origin to hit
                                                //self.set_depth_pixel(x as i32, y as i32, 0.0);
                                                break;
                                            }
                                        }
                                        None => {}
                                    }
                                }
                                None => {
                                    break;
                                }
                            }

                            let step_mult = 1.0;

                            if side_dist.x < side_dist.y {
                                if side_dist.x < side_dist.z {
                                    side_dist.x += delta_dist.x;
                                    map_pos.x += ray_step.x * step_mult;
                                    mask = BVec3::new(true, false, false);
                                } else {
                                    side_dist.z += delta_dist.z;
                                    map_pos.z += ray_step.z * step_mult;
                                    mask = BVec3::new(false, false, true);
                                }
                            } else {
                                if side_dist.y < side_dist.z {
                                    side_dist.y += delta_dist.y;
                                    map_pos.y += ray_step.y * step_mult;
                                    mask = BVec3::new(false, true, false);
                                } else {
                                    side_dist.z += delta_dist.z;
                                    map_pos.z += ray_step.z * step_mult;
                                    mask = BVec3::new(false, false, true);
                                }
                            }
                        }
                    }
                    _ => {}
                }
            }
        }
    }
}

//todo move edge things

pub fn edge_function<T>(a: &(T, T), b: &(T, T), c: &(T, T)) -> T
where
    T: Copy + Sub<Output = T> + Mul<Output = T>,
{
    (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
}

struct EdgeEquation {
    a: f32,
    b: f32,
    c: f32,
    tie: bool,
}

impl EdgeEquation {
    fn new(v0: &Vec3, v1: &Vec3) -> Self {
        let mut s = Self {
            a: v0.y - v1.y,
            b: v1.x - v0.x,
            c: 0.0,
            tie: false,
        };
        s.c = -(s.a * (v0.x + v1.x) + s.b * (v0.y + v1.y)) / 2.0;
        s.tie = if s.a != 0.0 { s.a > 0.0 } else { s.b > 0.0 };
        s
    }
    fn evaluate(&self, x: f32, y: f32) -> f32 {
        self.a * x + self.b * y + self.c
    }
    fn test_edge(&self, x: f32, y: f32) -> bool {
        self.test(self.evaluate(x, y))
    }
    fn test(&self, v: f32) -> bool {
        v > 0.0 || v == 0.0 && self.tie
    }
}

struct ParameterEquation {
    a: f32,
    b: f32,
    c: f32,
}

impl ParameterEquation {
    fn new(
        p0: f32,
        p1: f32,
        p2: f32,
        e0: &EdgeEquation,
        e1: &EdgeEquation,
        e2: &EdgeEquation,
        area: f32,
    ) -> Self {
        let factor = 1.0 / (2.0 * area);
        Self {
            a: factor * (p0 * e0.a + p1 * e1.a + p2 * e2.a),
            b: factor * (p0 * e0.b + p1 * e1.b + p2 * e2.b),
            c: factor * (p0 * e0.c + p1 * e1.c + p2 * e2.c),
        }
    }
    fn evaluate(&self, x: f32, y: f32) -> f32 {
        self.a * x + self.b * y + self.c
    }
}

//todo: rework line algarithm and move it somewere else
//todo: use dda to learn how to trace a grid
//todo:learn how to fo it in 3d
pub fn draw_line_pl(xa: i32, ya: i32, xb: i32, yb: i32) -> (Vec<(i32, i32)>, u8) {
    let mut flip = 0;
    let (mut x1, mut x2, mut y1, mut y2) = (xa, xb, ya, yb);
    let mut dx;
    let mut dy;
    let mut x;
    let mut y;
    let mut p;
    let mut point_list: Vec<(i32, i32)> = Vec::new();

    let m = if (x2 - x1) == 0 {
        y2 - y1
    } else {
        (y2 - y1) / (x2 - x1)
    };

    if m.abs() < 1 {
        if x1 > x2 {
            flip = 1;
            std::mem::swap(&mut x1, &mut x2);
            std::mem::swap(&mut y1, &mut y2);
        }

        dy = (y2 - y1).abs();
        dx = (x2 - x1).abs();

        p = 2 * dy - dx;

        y = y1;
        x = x1;

        while x <= x2 {
            point_list.push((x, y));

            x = x + 1;
            if p >= 0 {
                if m >= 1 || y1 < y2 {
                    y = y + 1;
                } else {
                    y = y - 1;
                }
                p = p + 2 * dy - 2 * dx;
            } else {
                p += 2 * dy;
            }
        }
    }
    if m.abs() >= 1 {
        {
            if y1 > y2 {
                flip = 1;
                std::mem::swap(&mut x1, &mut x2);
                std::mem::swap(&mut y1, &mut y2);
            }

            dy = (y2 - y1).abs();
            dx = (x2 - x1).abs();

            p = 2 * dx - dy;

            y = y1;
            x = x1;

            while y < y2 {
                point_list.push((x, y));
                y += 1;
                p = if p >= 0 {
                    x = if m >= 1 { x + 1 } else { x - 1 };
                    p + 2 * dx - 2 * dy
                } else {
                    p + 2 * dx
                };
            }
        }
    }
    (point_list, flip)
}

/// used to help pass on custom values between stages
pub struct SharderPipe<T, V, F>
where
    V: Fn(&mut Vertex, &mut T),
    F: Fn(Vertex, &mut T) -> Pixel,
{
    pub back_side_cull: bool,
    pub data: T,
    vertex_shader: V,
    frag_shader: F,
}

impl<T, V, F> SharderPipe<T, V, F>
where
    V: Fn(&mut Vertex, &mut T),
    F: Fn(Vertex, &mut T) -> Pixel,
{
    pub fn new(data: T, vertex_shader: V, frag_shader: F) -> Self {
        Self {
            back_side_cull: true,
            data,
            vertex_shader,
            frag_shader,
        }
    }

    fn run_vertex_shader(&mut self, vert: &mut Vertex) {
        (self.vertex_shader)(vert, &mut self.data);
    }
    fn run_frag_shader(&mut self, vert: Vertex) -> Pixel {
        (self.frag_shader)(vert, &mut self.data)
    }
}
