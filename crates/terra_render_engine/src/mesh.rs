extern crate gltf;
extern crate prisma_maths;

use rayon::prelude::*;

use dot_vox::*;
//use gltf::Gltf;
use prisma_maths::*;

pub enum MeshRenderMode {
    Dots,
    VoxelPoint,
    Lines,
    WireFrame,
    Tri,
}

#[derive(Clone, Copy)]
pub struct Vertex {
    pub point: Vec4,
    pub normal: Vec3,
    pub color: Vec3,
    pub uv: Vec2,
}

impl Vertex {
    pub fn new(point: Vec4, normal: Vec3, color: Vec3, uv: Vec2) -> Self {
        Self {
            point,
            normal,
            color,
            uv,
        }
    }
    pub fn leap(&self, other: Vertex, s: f32) -> Vertex {
        Vertex::new(
            self.point.lerp(other.point, s),
            self.normal.lerp(other.normal, s),
            self.color.lerp(other.color, s),
            self.uv.lerp(other.uv, s),
        )
    }
}

#[derive(Clone, Copy)]
pub struct Triangle {
    pub p: (Vertex, Vertex, Vertex),
}

impl Triangle {
    pub fn new(p: (Vertex, Vertex, Vertex)) -> Self {
        Self { p }
    }
}

pub struct Material {
    // vertex shader
// pixel/frag shader
}
pub struct Mesh {
    //name: String,
    pub transform: Mat4,
    pub vertices: Vec<Vertex>,
    pub v_index: Vec<usize>, // uses to index vertices, every 3 v_index is a tri
}

impl Mesh {
    pub fn new() -> Mesh {
        Mesh {
            transform: Mat4::identity(),
            vertices: Vec::new(),
            v_index: Vec::new(),
        }
    }
    pub fn load_glft_mesh(path: &String) -> Mesh {
        let mut o_mesh = Mesh::new();
        let (gltf, buffers, _) = gltf::import(path).unwrap();

        for mesh in gltf.meshes() {
            //println!("Mesh #{}", mesh.index());

            for primitive in mesh.primitives() {
                //println!("- Primitive #{}", primitive.index());
                //let idx = primitive.indices().unwrap();

                let reader = primitive.reader(|buffer| Some(&buffers[buffer.index()]));

                if let Some(iter) = reader.read_indices() {
                    for i in iter.into_u32() {
                        o_mesh.v_index.push(i as usize);
                    }
                }

                let mut tpos = Vec::new();
                if let Some(iter) = reader.read_positions() {
                    for vertex_position in iter {
                        //println!("{:?}", vertex_position);
                        let v = Vec4::new(
                            vertex_position[0],
                            vertex_position[1],
                            vertex_position[2],
                            1.,
                        );

                        tpos.push(v);
                    }
                }
                let mut tnorm = Vec::new();
                if let Some(iter) = reader.read_normals() {
                    for vertex_position in iter {
                        //println!("{:?}", vertex_position);
                        let v =
                            Vec3::new(vertex_position[0], vertex_position[1], vertex_position[2]);

                        tnorm.push(v);
                    }
                }
                let mut tcol = Vec::new();
                if let Some(iter) = reader.read_colors(0) {
                    for vertex_position in iter.into_rgb_f32() {
                        //println!("{:?}", vertex_position);
                        let v =
                            Vec3::new(vertex_position[0], vertex_position[1], vertex_position[2]);

                        tcol.push(v);
                    }
                }
                let mut tuv = Vec::new();
                if let Some(iter) = reader.read_tex_coords(0) {
                    for vertex_position in iter.into_f32() {
                        //println!("{:?}", vertex_position);
                        let v = Vec2::new(vertex_position[0], vertex_position[1]);

                        tuv.push(v);
                    }
                }
                for i in 0..tpos.len() {
                    o_mesh
                        .vertices
                        .push(Vertex::new(tpos[i], tnorm[i], tcol[i], tuv[i]));
                }
            }
        }

        o_mesh
    }
}

pub struct VoxelTexture {
    pub map: Vec<u8>,
    pub distance_field: Vec<u8>,
    pub palette: Vec<u32>,
    pub size: (u32, u32, u32),
}

impl VoxelTexture {
    pub fn new(map: Vec<u8>, size: (u32, u32, u32)) -> Self {
        let ds = vec![255; map.len()];

        Self {
            map,
            distance_field: ds,
            size,
            palette: Vec::new(),
        }
    }
    pub fn create_voxel_texture(size: (u32, u32, u32)) -> VoxelTexture {
        let mut map = Vec::new();
        for i in 0..(size.0 * size.1 * size.2) {
            map.push(0);
        }
        VoxelTexture::new(map, size)
    }
    pub fn get_voxel_index(&self, x: i32, y: i32, z: i32) -> Option<i32> {
        if x < self.size.0 as i32
            && x >= 0
            && y < self.size.1 as i32
            && y >= 0
            && z < self.size.2 as i32
            && z >= 0
        {
            return Some(
                (z * self.size.0 as i32 * self.size.1 as i32) + (y * self.size.0 as i32) + x,
            );
        } else {
            return None;
        }
    }
    /*
    pub fn get_index_voxel(&self, index: usize) -> (u8, u8, u8) {
        let mut idx = index as u8;
        let z = idx / (self.size.0 * self.size.1);
        idx -= z * self.size.0 * self.size.1;
        let y = idx / self.size.0;
        let x = idx % self.size.0;
        (x, y, z)
    }
    */
    pub fn set_voxel(&mut self, x: i32, y: i32, z: i32, idx: u8) {
        let i = self.get_voxel_index(x, y, z);
        match i {
            Some(id) => {
                self.map[id as usize] = idx;
            }
            _ => {}
        }
    }
    pub fn get_voxel(&mut self, x: i32, y: i32, z: i32) -> Option<u8> {
        let i = self.get_voxel_index(x, y, z);
        match i {
            Some(id) => Some(self.map[id as usize]),
            _ => None,
        }
    }
    pub fn clear_volume(&mut self) {
        for v in self.map.iter_mut() {
            *v = 0;
        }
    }
    pub fn tm(&mut self, other: &mut VoxelTexture, tm: &Mat4) {
        let v1 = *tm * Vec4::new(0.0, 0.0, 0.0, 1.0);
        let v2 = *tm * Vec4::new(other.size.0 as f32, 0.0, 0.0, 1.0);
        let v3 = *tm * Vec4::new(0.0, other.size.1 as f32, 0.0, 1.0);
        let v4 = *tm * Vec4::new(other.size.0 as f32, other.size.1 as f32, 0.0, 1.0);

        let v5 = *tm * Vec4::new(0.0, 0.0, other.size.2 as f32, 1.0);
        let v6 = *tm * Vec4::new(other.size.0 as f32, 0.0, other.size.2 as f32, 1.0);
        let v7 = *tm * Vec4::new(0.0, other.size.1 as f32, other.size.2 as f32, 1.0);
        let v8 = *tm
            * Vec4::new(
                other.size.0 as f32,
                other.size.1 as f32,
                other.size.2 as f32,
                1.0,
            );

        let min_v = v1.min(v2.min(v3.min(v4.min(v5.min(v6.min(v7.min(v8)))))));
        let max_v = v1.max(v2.max(v3.max(v4.max(v5.max(v6.max(v7.max(v8)))))));

        let inverse_m = tm.inverse(); // need to do loop in side ove sprite valume to render volume

        for z in (min_v.z as i32..max_v.z as i32).rev() {
            for x in (min_v.x as i32..max_v.x as i32).rev() {
                for y in (min_v.y as i32..max_v.y as i32).rev() {
                    for offx in 0..3 {
                        for offy in 0..3 {
                            let offval = 0.7;
                            let off_set_x = offx as f32 * offval;
                            let off_set_y = offy as f32 * offval;
                            // sample form pos
                            let n_v = (inverse_m
                                * Vec4::new(
                                    x as f32 + off_set_x,
                                    y as f32 + off_set_y,
                                    z as f32 + off_set_x,
                                    1.0,
                                ))
                            .round();

                            match other.get_voxel(n_v.x as i32, n_v.y as i32, n_v.z as i32) {
                                Some(v) => {
                                    // sameple to pos
                                    let mut tfv = Vec4::new(x as f32, y as f32, z as f32, 1.0);
                                    tfv = tfv.round();

                                    self.set_voxel(
                                        (tfv.x + off_set_x) as i32,
                                        (tfv.y + off_set_y) as i32,
                                        (tfv.z + off_set_x) as i32,
                                        v,
                                    );

                                    // may have to floor
                                }
                                None => {}
                            }
                        }
                    }
                }
            }
        }
    }
    pub fn load_voxel_map_from_vox(path: &String) -> VoxelTexture {
        let vm = load(path).unwrap();
        let mut voxel_map = VoxelTexture::new(Vec::new(), (0, 0, 0));
        voxel_map.palette = vm.palette;
        for model in vm.models.iter() {
            voxel_map.size = (model.size.x, model.size.z, model.size.y);
            for _ in 0..(voxel_map.size.0 as usize)
                * (voxel_map.size.1 as usize)
                * (voxel_map.size.2 as usize)
            {
                voxel_map.map.push(0);
            }

            for vox in model.voxels.iter() {
                let idx = voxel_map.get_voxel_index(vox.x as i32, vox.z as i32, vox.y as i32);
                match idx {
                    Some(id) => {
                        voxel_map.map[id as usize] = vox.i;
                    }
                    _ => {}
                }
            }
        }
        voxel_map.update_distance_filed();
        voxel_map
    }

    pub fn update_distance_filed(&mut self) {
        let ds = vec![8; self.map.len()];
        self.distance_field = ds;
        for x in 0..self.size.0 {
            //println!("{:?}", x);
            for y in 0..self.size.1 {
                for z in 0..self.size.2 {
                    let vp = (x, y, z);

                    let mut dist = 255;
                    let current_v = self.get_voxel(x as i32, y as i32, z as i32).unwrap();
                    let current_v_index =
                        self.get_voxel_index(x as i32, y as i32, z as i32).unwrap() as usize;

                    if current_v != 0 {
                        self.distance_field[current_v_index] = 0;
                        continue;
                    };

                    let range = 8;
                    for xi in
                        ((-range) + x as i32).max(0)..(range + x as i32).min(self.size.0 as i32)
                    {
                        for yi in
                            ((-range) + y as i32).max(0)..(range + y as i32).min(self.size.1 as i32)
                        {
                            for zi in ((-range) + z as i32).max(0)
                                ..(range + z as i32).min(self.size.2 as i32)
                            {
                                if dist == 0 {
                                    self.distance_field[current_v_index] = 0;
                                    continue;
                                }
                                let n_v = self.get_voxel(xi as i32, yi as i32, zi as i32).unwrap();

                                if n_v != 0 {
                                    let new_vp = (xi, yi, zi);

                                    let n_dist = vec3(vp.0 as f32, vp.1 as f32, vp.2 as f32)
                                        .distance_squared(vec3(
                                            new_vp.0 as f32,
                                            new_vp.1 as f32,
                                            new_vp.2 as f32,
                                        )) as u8;
                                    if n_dist < dist {
                                        dist = n_dist;
                                        self.distance_field[current_v_index] = dist;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
