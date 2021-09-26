extern crate glam as glm;

// expand glm to have affine transformations, geting a afine matrix or talking a matrix and apllying a afine to it

mod math_functions;
mod types;

pub use glam::*;
pub use math_functions::*;
pub use types::*;

// need to be abule to offset pos relative to rot ie fps,
pub struct CamTransform {}

pub struct btransform {
    pos: Vec3,
    rot: Quat,
    scale: Vec3,
}

pub fn perspective_proj(fov: f32, wight: f32, hight: f32, z_near: f32, z_far: f32) -> Mat4 {
    //a.x = (a.x / (ar * a.z * tanhfov)) * hw + hw;
    //        a.y = (a.y / (a.z * tanhfov)) * hh + hh;

    let ar = wight / hight;
    let t = ((fov / 2.).to_radians()).tan();
    let z_range = z_near - z_far;

    mat4(
        vec4(1. / ar * t, 0.0, 0.0, 0.0),
        vec4(0.0, 1. / t, 0.0, 0.0),
        vec4(
            0.0,
            0.0,
            -z_near - z_far / z_range,
            2.0 * z_far * z_near / z_range,
        ),
        vec4(0.0, 0.0, 1.0, 0.0),
    )
}

pub fn translate_mat3(mp: &mut glm::Mat3, x: f32, y: f32) {
    mp.z_axis.x += x;
    mp.z_axis.y += y;
}

//bugd
pub fn rotate_mat3(mp: &mut glm::Mat3, x: f32) {
    mp.x_axis.x = x.cos();
    mp.x_axis.y = -x.sin();
    mp.y_axis.x = x.sin();
    mp.y_axis.y = x.cos();
}

pub fn scale_mat3(mp: &mut glm::Mat3, w: f32, h: f32) {
    mp.x_axis.x += w;
    mp.y_axis.y += h;
}

pub fn shear_mat3(mp: &mut glm::Mat3, x: f32, y: f32) {
    mp.y_axis.x += x;
    mp.x_axis.y += y;
}

#[cfg(test)]
mod tests {
    use super::*;
    //use num_traits::Pow;
    use types::Vector3;

    
    #[test]
    fn remap_test_i32() {
        assert_eq!(remap(50, 0, 100, 200, 400), 300);
    }
    #[test]
    fn remap_test_f32() {
        assert_eq!(remap(0.5, 0., 1., 0., 1000.), 500.);
    }
    #[test]
    fn remapi32dp_test() {
        assert_eq!(remap_i32dp(255, 0, 255, 255, 0), 0);
    }

    #[test]
    fn vector3_dot() {
        let a = Vector3::new(5, 5, 5);
        let b = Vector3::new(-5, 5, -0);
        let c = Vector3::new(7, 3, 0);
        assert_eq!(a.dot(&b), 0);
        assert_eq!(a.dot(&c), 50);
    }

    #[test]
    fn vector3_cross() {
        let a = Vector3::new(-1, -2, 3);
        let b = Vector3::new(4, 0, -8);
        let c = Vector3::new(7, 3, 1);
        assert_eq!(a.cross(&b), Vector3::new(16, 4, 8));
        assert_eq!(a.cross(&c), Vector3::new(-11, 22, 11));
    }
    #[test]
    fn vector3_length() {
        let _a = Vector3::new(-1, -2, 3);
        let c = Vector3::new(5.0, 5.0, 5.0);

        //assert_eq!(a.length(), 3);
        assert_eq!(c.length() as i32, 8.66025 as i32);
    }

    #[test]
    fn vector3_normalize() {
        let _a = Vector3::new(1., 1., 0.0);
        let _b = Vector3::new(0.70710677, 0.70710677, 0.0);
        let _c = Vector3::new(7.0, 3.0, 1.0);
        //assert_eq!(c.length(), 11.0);
    }
    #[test]
    fn vector3_abs() {
        let a = Vector3::new(1., -0.6, 0.0);
        let b = Vector3::new(4, 0, -8);
        let c = Vector3::new(-7.0, -3.0, 1.0);
        assert_eq!(a.abs(), Vector3::new(1., 0.6, 0.0));
        assert_eq!(b.abs(), Vector3::new(4, 0, 8));
        assert_eq!(c.abs(), Vector3::new(7.0, 3.0, 1.0));
    }
}
