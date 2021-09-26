extern crate prisma_maths;

use super::ray_tracing::*;
use prisma_maths::*;

//todo: inverse_camera bools
#[derive(Clone, Copy)]
pub struct Camera {
    pub position: Vec3,
    pub rotation: Quat,
    pub fov: f32,
    pub aspect_ratio: f32,
    pub z_near: f32,
    pub z_far: f32,
    dirty_camera_to_world: bool, // pos or rot update?
    dirty_perspective: bool,     //
    camera_to_world_mat: Mat4,
    perspective_mat: Mat4,
    inverse_camera_to_world_mat: Mat4,
    inverse_perspective_mat: Mat4,
}

impl Camera {
    pub fn new(fov: f32, aspect_ratio: f32, z_near: f32, z_far: f32) -> Self {
        Self {
            position: Vec3::new(0., 0., 0.),
            rotation: Quat::from_rotation_z(0.0),
            fov,
            aspect_ratio,
            z_near,
            z_far,
            dirty_camera_to_world: true,
            dirty_perspective: true,
            perspective_mat: Mat4::identity(),
            camera_to_world_mat: Mat4::identity(),
            inverse_perspective_mat: Mat4::identity(),
            inverse_camera_to_world_mat: Mat4::identity(),
        }
    }
    pub fn set_position(&mut self, pos: Vec3) {
        self.position = pos;
        self.dirty_camera_to_world = true;
    }
    pub fn set_rotation(&mut self, rot: Quat) {
        self.rotation = rot;
        self.dirty_camera_to_world = true;
    }
    pub fn get_camera_to_world(&mut self) -> Mat4 {
        if self.dirty_camera_to_world {
            self.dirty_camera_to_world = false;
            self.camera_to_world_mat =
                Mat4::from_quat(self.rotation) * Mat4::from_translation(self.position);
            self.inverse_camera_to_world_mat = self.camera_to_world_mat.inverse();
        }
        self.camera_to_world_mat
    }

    pub fn get_inverse_camera_to_world(&mut self) -> Mat4 {
        self.inverse_camera_to_world_mat
    }

    pub fn get_perspective_mat4(&mut self) -> Mat4 {
        if self.dirty_perspective {
            self.dirty_perspective = false;
            self.perspective_mat = Mat4::perspective_rh_gl(
                self.fov.to_radians(),
                self.aspect_ratio,
                self.z_near,
                self.z_far,
            );
            //self.perspective_mat = Mat4::orthographic_rh_gl(-1.0, 1.0, -1.0, 1.0, 1.0, 100.0);
            self.inverse_perspective_mat = self.perspective_mat.inverse();
        }
        self.perspective_mat
    }

    pub fn get_inverse_perspective_mat4(&mut self) -> Mat4 {
        self.inverse_perspective_mat
    }

    pub fn get_camera_ray(&mut self, uv: Vec2) -> Ray {
        let origin = (self.get_inverse_camera_to_world() * vec4(0., 0., 0., 1.)).xyz();

        let mut direction = (self.get_inverse_perspective_mat4() * vec4(uv.x, uv.y, 0., 1.)).xyz();

        direction = (self.get_inverse_camera_to_world()
            * Vec4::new(direction.x, direction.y, direction.z, 0.))
        .xyz();

        direction = direction.normalize();

        Ray::create_ray(origin, direction)
    }
}
