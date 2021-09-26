extern crate prisma_maths;

use prisma_maths::*;

use crate::VoxelTexture;

pub struct Ray {
    pub origin: Vec3,
    pub direction: Vec3,
    pub energy: Vec3,
}

pub struct RayHit {
    pub position: Vec3,
    pub distance: f32,
    pub normal: Vec3,
}

impl RayHit {
    pub fn new(position: Vec3, distance: f32, normal: Vec3) -> Self {
        Self {
            position,
            distance,
            normal,
        }
    }

    pub fn create_ray_hit() -> RayHit {
        return RayHit::new(
            Vec3::new(0.0, 0.0, 0.0),
            f32::INFINITY,
            Vec3::new(0.0, 0.0, 0.0),
        );
    }
}

impl Ray {
    pub fn new(origin: Vec3, direction: Vec3, energy: Vec3) -> Self {
        Self {
            origin,
            direction,
            energy,
        }
    }
    pub fn create_ray(origin: Vec3, direction: Vec3) -> Self {
        Self {
            origin,
            direction,
            energy: Vec3::new(1., 1., 1.),
        }
    }
}

const STEPS: i32 = 256;
//const GRID: f32 = 16.0;
const EPSY: f32 = 0.001;

// ray navagating a voxel grid
// ref: https://www.shadertoy.com/view/ltSXDm
pub fn ray_voxel_trace(ray: &Ray, voxel_grid: &VoxelTexture) -> Option<(f32, f32, f32)> {
    let origin = ray.origin;
    let direction = ray.direction;
    let mut point = origin.floor();
    let steps = Vec3::new(sign(direction.x), sign(direction.y), sign(direction.z));
    let t_delta = 1.0 / direction.abs().max(Vec3::splat(EPSY));
    let mut t_max: Vec3 = (Vec3::splat(1.0)
        - Vec3::new(
            (origin.x * steps.x).fract(),
            (origin.y * steps.y).fract(),
            (origin.z * steps.z).fract(),
        ))
        * t_delta;

    for _ in 0..STEPS {
        // voxel_look up (point);
        let voxid = voxel_grid.get_voxel_index(point.x as i32, point.y as i32, point.z as i32);

        match voxid {
            Some(vi) => {
                let vox = *voxel_grid.map.get(vi as usize).unwrap();
                if vox != 0 {
                    return Some((point.x as f32, point.y as f32, point.z as f32));
                //break;
                } else {
                    //nothing
                }
            }
            None => {}
        }

        let a = t_max.cmplt(t_max.yzx()); // may have to rework
        let b = t_max.cmple(t_max.zxy());
        let select: BVec3 = a & b; // may not work

        point += Vec3::select(select, steps, Vec3::splat(0.0));
        t_max += Vec3::select(select, t_delta, Vec3::splat(0.0));
    }
    None
}

pub fn interest_sphere(ray: &Ray, best_hit: &mut RayHit, sphere: Vec4) {
    let d = ray.origin - sphere.xyz();
    let p1 = -ray.direction.dot(d);
    let p2sqr = p1 * p1 - d.dot(d) + sphere.w * sphere.w;
    if p2sqr < 0.0 {
        return;
    }
    let p2 = p2sqr.sqrt();
    let t = if p1 - p2 > 0.0 { p1 - p2 } else { p1 + p2 };
    if t > 0.0 && t < best_hit.distance {
        best_hit.distance = t;
        best_hit.position = ray.origin + t * ray.direction;
        best_hit.normal = (best_hit.position - sphere.xyz()).normalize();
    }
}

/// Calcs intersection and exit distances, normal, face and UVs,
/// row is the ray origin in world space,
/// rdw is the ray direction in world space,
/// txx is the world-to-box transformation,
/// txi is the box-to-world transformation,
/// ro and rd are in world space,
/// rad is the half-length of the box,
///
/// oT contains the entry and exit points,
/// oN is the normal in world space,
/// oU contains the UVs at the intersection point,
/// oF contains the index if the intersected face [0..5],
pub fn box_intersect(
    row: &Vec3,
    rdw: &Vec3,
    txx: &Mat4,
    txi: &Mat4,
    rad: &Vec3,
    ot: &mut Vec2,
    on: &mut Vec3,
    ou: &mut Vec2,
    of: &mut i32,
) -> bool {
    // convert from world to box space
    let rd = (*txx * rdw.extend(0.0)).xyz();
    let ro = (*txx * row.extend(1.0)).xyz();

    // ray-box intersection in box space
    let m = 1.0 / rd;
    let s = Vec3::new(
        if rd.x < 0.0 { 1.0 } else { -1.0 },
        if rd.y < 0.0 { 1.0 } else { -1.0 },
        if rd.z < 0.0 { 1.0 } else { -1.0 },
    );
    let t1 = m * (-ro + s * *rad);
    let t2 = m * (-ro - s * *rad);
    let tn = max!(t1.x, t1.y, t1.z);
    let tf = min!(t2.x, t2.y, t2.z);

    if tn > tf || tf < 0.0 {
        return false;
    }

    // compute normal (in world space), face and UV
    if t1.x > t1.y && t1.x > t1.z {
        *on = txi.x_axis.xyz() * s.x;
        *ou = ro.yz() + rd.yz() * t1.x;
        *of = (1 + s.x as i32) / 2;
    } else if t1.y > t1.z {
        *on = txi.y_axis.xyz() * s.y;
        *ou = ro.zx() + rd.zx() * t1.y;
        *of = (5 + s.y as i32) / 2;
    } else {
        *on = txi.z_axis.xyz() * s.z;
        *ou = ro.xy() + rd.xy() * t1.z;
        *of = (9 + s.z as i32) / 2;
    }
    *ot = Vec2::new(tn, tf);

    return true;
}
