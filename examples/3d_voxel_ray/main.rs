use pixels::{Pixels, SurfaceTexture};
use winit::window::WindowBuilder;
use winit::{
    dpi::PhysicalPosition,
    event_loop::{ControlFlow, EventLoop},
};
use winit::{
    dpi::Position,
    event::{Event, VirtualKeyCode, WindowEvent},
};
use winit_input_helper::WinitInputHelper;

use prisma_maths::*;
use terra_render_engine as tre;
use tre::{Camera, SharderPipe, VoxelTexture};
use windowing::*;

const BUFF_SCALE: u32 = 18;
const BUFF_WIDTH: u32 = 16 * BUFF_SCALE;
const BUFF_HEIGHT: u32 = 9 * BUFF_SCALE;

fn main() {
    let mut r = tre::RenderEngine::new((BUFF_WIDTH as usize, BUFF_HEIGHT as usize));
    r.init_depth_buffer();

    // load texture and get size

    let font_bit = tre::RetroBitMap::load_bitmap_from_image("assets/textures/ass_font_tran.png");

    // select what vox model to load
    let voxel_texture = match 10 {
        0 => VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/test box map.vox")),
        1 => VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/chr_fox.vox")),
        2 => VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/chr_old.vox")),
        3 => VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/chr_sword.vox")),
        4 => VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/chr_knight.vox")),
        5 => VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/doom.vox")),
        6 => VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/monu1.vox")),
        7 => VoxelTexture::load_voxel_map_from_vox(&String::from(
            "assets/vox/monu6-without-water.vox",
        )),
        8 => VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/monu8.vox")),
        9 => VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/monu9.vox")),
        10 => VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/monu10.vox")),
        11 => VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/monu16.vox")),
        12 => VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/menger.vox")),
        13 => VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/nature.vox")),
        14 => VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/room.vox")),
        15 => {
            VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/red_booth_solid.vox"))
        }
        16 => VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/street_scene.vox")),
        17 => VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/treehouse.vox")),
        _ => panic!("o no"),
    };
    let voxe_model2 =
        VoxelTexture::load_voxel_map_from_vox(&String::from("assets/vox/chr_fox.vox"));

    //let mut voxel_render_volume = VoxelTexture::create_voxel_texture((128, 128, 128));
    //voxel_render_volume.palette = voxe_model2.palette.clone();

    // select what mesh model to load
    let mut mp = Mat4::identity();
    mp.w_axis.z = 0.0;
    let mut model = match 5 {
        0 => {
            mp.w_axis.z = -2.0;
            tre::Mesh::load_glft_mesh(&String::from("assets/models/suzanne.glb"))
        }
        1 => {
            mp.w_axis.z = -5.0;
            tre::Mesh::load_glft_mesh(&String::from("assets/models/plane.glb"))
        }
        2 => {
            mp.w_axis.z = -4.0;
            tre::Mesh::load_glft_mesh(&String::from("assets/models/box_tri.glb"))
        }
        3 => {
            // 20_000 tri sphere
            mp.w_axis.z = -1.5;
            tre::Mesh::load_glft_mesh(&String::from("assets/models/hp_ico_sphere.glb"))
        }
        4 => {
            // 122_000 tri sphere
            mp.w_axis.z = -1.5;
            tre::Mesh::load_glft_mesh(&String::from("assets/models/shp_ico.glb"))
        }
        5 => {
            mp.w_axis.z = -10.0;
            tre::Mesh::load_glft_mesh(&String::from("assets/models/suzanne_ss.glb"))
        }
        6 => {
            mp.w_axis.z = 5.0;
            tre::Mesh::load_glft_mesh(&String::from("assets/models/tri.glb"))
        }
        _ => tre::Mesh::new(),
    };
    let _terrain_model = tre::Mesh::load_glft_mesh(&String::from("assets/models/terrain.glb"));

    let mut input = WinitInputHelper::new();
    let event_loop = EventLoop::new();
    let window = WindowBuilder::new()
        .with_title("rust by examp !")
        .build(&event_loop)
        .unwrap();
    //window.set_cursor_visible(false);
    //window.set_cursor_grab(true);

    let window_size = window.inner_size();
    let surface_texture = SurfaceTexture::new(window_size.width, window_size.height, &window);
    let mut pixels = Pixels::new(BUFF_WIDTH, BUFF_HEIGHT, surface_texture).unwrap();

    let mut camera = Camera::new(90.0, 16. / 9., 1.0, 40.0);
    let mut cam_rot = (0.0, 0.0);
    let mut cam_vel = vec3(0., 0., 0.);

    let hw = r.buff_size.0 as f32 / 2.;
    let hh = r.buff_size.1 as f32 / 2.;
    let mut event_timer = FrameTimer::new(0.04);
    let mut render_timer = FrameTimer::new(0.04);

    let mut m_pos = (0.0_f32, 0.0_f32);

    event_loop.run(move |event, _, control_flow| {
        ////////////
        input.update(&event);

        match event {
            Event::NewEvents(_) => {}
            Event::WindowEvent {
                window_id: _,
                event,
            } => match event {
                WindowEvent::Resized(s) => pixels.resize_surface(s.width, s.height),
                _ => {}
            },
            Event::DeviceEvent {
                device_id: _,
                event: _,
            } => {}
            Event::UserEvent(_) => {}
            Event::Suspended => {}
            Event::Resumed => {}
            Event::MainEventsCleared => {
                event_timer.start();
                let n_window_size = window.inner_size();
                let mouse_sens = 0.04 * render_timer.delta_time as f32;

                let mut m_dif = (0.0, 0.0);
                match input.mouse() {
                    Some(pos) => {
                        if (pos.0 * pos.1) != 0.0 {
                            m_dif.0 = m_pos.0 - pos.0;
                            m_dif.1 = m_pos.1 - pos.1;

                            m_pos = (
                                (n_window_size.width as f32 / 2.).floor(),
                                (n_window_size.height as f32 / 2.).floor(),
                            );
                            window
                                .set_cursor_position(Position::Physical(PhysicalPosition::new(
                                    n_window_size.width as i32 / 2,
                                    n_window_size.height as i32 / 2,
                                )))
                                .unwrap();
                        }
                    }
                    None => {}
                }

                cam_rot.0 += m_dif.0 * mouse_sens;
                cam_rot.1 += m_dif.1 * mouse_sens;
                camera.set_rotation(
                    Quat::from_rotation_ypr(0.0, -cam_rot.1, 0.0)
                        * Quat::from_rotation_ypr(-cam_rot.0, 0.0, 0.0),
                );
                let mut imput_v = vec3(0., 0., 0.);

                if input.key_held(VirtualKeyCode::A) {
                    imput_v.x += 1.0;
                }
                if input.key_held(VirtualKeyCode::D) {
                    imput_v.x -= 1.0;
                }
                if input.key_held(VirtualKeyCode::W) {
                    imput_v.z += 1.0;
                }
                if input.key_held(VirtualKeyCode::S) {
                    imput_v.z -= 1.0;
                }

                imput_v = Mat4::from_quat(camera.rotation)
                    .inverse()
                    .transform_point3(imput_v);

                if input.key_held(VirtualKeyCode::Q) {
                    imput_v.y += 1.0;
                }
                if input.key_held(VirtualKeyCode::E) {
                    imput_v.y -= 1.0;
                }

                let frc = 0.0465;
                // velocity -= min(abs(velocity.length()), frc) * sign(velocity.length()) * velocity.normalized() * delta
                if cam_vel.length() != 0.0 {
                    cam_vel -= (cam_vel.length().abs().min(frc) * sign(cam_vel.length()))
                        * cam_vel.normalize();
                }
                let max_speed = if !input.key_held(VirtualKeyCode::LShift) {
                    if !input.key_held(VirtualKeyCode::LControl) {
                        2.0
                    } else {
                        0.4
                    }
                } else {
                    25.0
                };

                cam_vel += imput_v * 1.0;
                if cam_vel.length() > max_speed {
                    cam_vel = cam_vel.normalize() * (max_speed - 0.0001);
                }

                camera.set_position(camera.position + (cam_vel * render_timer.delta_time as f32));
                if input.key_released(VirtualKeyCode::Escape) {
                    *control_flow = ControlFlow::Exit;
                    return;
                }
                window.request_redraw();
                event_timer.end();
            }
            Event::RedrawRequested(_) => {
                render_timer.start();
                r.clear_depth_buffer();

                let cdata = (0.0_f32, 0);
                let mut standerd_mat = SharderPipe::new(
                    cdata,
                    |_v, _d| {},
                    |v, _d| {
                        let light_level = clamp(
                            v.normal.dot(vec3(1.0, 1.0, 1.0).normalize()) * 0.6 + 0.4,
                            0.0,
                            1.0,
                        ) * 2.0;
                        // todo add fn to pixel to return a float tupel as pixel
                        // R = max(0, min(1, bias + scale * (1.0 + I • N)power))

                        tre::Pixel::new(
                            (v.color.x * light_level * 255.0) as u8,
                            (v.color.y * light_level * 255.0) as u8,
                            (v.color.z * light_level * 255.0) as u8,
                            255,
                        )
                    },
                );

                for i in 0..BUFF_WIDTH * BUFF_HEIGHT {
                    let (x, y) = tre::index_to_xy(i, BUFF_WIDTH);
                    let u = remap(x, 0, BUFF_WIDTH, 0, 100) as u8;
                    let v = remap(y, 0, BUFF_HEIGHT, 0, 100) as u8;

                    r.buffer
                        .set_pixel(x as i32, y as i32, &tre::Pixel::new(u, v, 80, 255));
                }

                let pos = Vec4::new(
                    (event_timer.total_time as f32 * 0.7).sin() * 8.,
                    (((event_timer.total_time as f32 * 0.6).sin() * 2.3).cos() * 4.) + 2.,
                    ((event_timer.total_time as f32).cos() * 8.) - 15.,
                    1.,
                );

                mp.w_axis = pos;
                model.transform = mp
                    * Mat4::from_axis_angle(
                        Vec3::new(2., 5., 1.).normalize(),
                        event_timer.total_time as f32,
                    );
                r.draw_mesh(
                    &model,
                    &mut standerd_mat,
                    &mut camera,
                    &tre::MeshRenderMode::Tri,
                );

                mp.w_axis = Vec4::new(1., 1., 10., 1.);
                model.transform =
                    mp * Mat4::from_axis_angle(Vec3::new(0.0, 0.0, 1.).normalize(), 0.0);
                r.draw_mesh(
                    &model,
                    &mut standerd_mat,
                    &mut camera,
                    &tre::MeshRenderMode::Tri,
                );

                for i in 0..3 {
                    for j in 0..3 {
                        let pos = Vec4::new(
                            i as f32 * 5. + 5.,
                            ((((event_timer.total_time as f32 * 0.2) + i as f32 + j as f32).sin()
                                * 1.1431)
                                .cos()
                                + 1.0)
                                * 4.,
                            j as f32 * -5.,
                            1.,
                        );
                        mp.w_axis = pos;
                        model.transform = mp
                            * Mat4::from_axis_angle(
                                Vec3::new(2., 5., 1.).normalize(),
                                event_timer.total_time as f32
                                    * ((i as f32 * 0.4) + (j as f32 * 0.8)),
                            );
                        //r.draw_mesh(&model, &mut camera, &tre::MeshRenderMode::Tri);
                    }
                }

                //r.draw_mesh(&terrain_model, &mut camera, &tre::MeshRenderMode::Tri);

                // ray trace
                if true {
                    let pm = Mat4::from_translation(Vec3::new(0., 0., 4.));
                    let tm = pm * Mat4::mul_scalar(&mut Mat4::identity(), 1.0);
                    //tm = Mat4::from_axis_angle(
                    //    Vec3::new(0., 1., 0.).normalize(),
                    //    (33.0 as f32).to_radians(),
                    //) * tm;

                    r.voxel_ray_trace(&voxel_texture, &tm, &mut camera);

                    let pm = Mat4::from_translation(Vec3::new(4., 0., -4.));
                    let mut tm = pm;
                    tm = Mat4::from_axis_angle(
                        Vec3::new(3.0, 1., 2.).normalize(),
                        render_timer.total_time as f32,
                    ) * tm;

                    r.voxel_ray_trace(&voxe_model2, &tm, &mut camera);
                }
                let cpos = (camera.get_camera_to_world().inverse() * vec4(0., 0., 0., 1.)).xyz();

                r.buffer.draw_string(
                    &font_bit,
                    (-0.95 * hw + hw) as i32,
                    (0.70 * hh + hh) as i32,
                    (8, 8),
                    &format!(
                        "CAM POS:x{:#0.01} y{:#0.01} z{:#0.01}",
                        cpos.x, cpos.y, cpos.z
                    ),
                    Some(&tre::Pixel::new(230, 75, 10, 200)),
                );

                r.buffer.draw_string(
                    &font_bit,
                    (-0.95 * hw + hw) as i32,
                    (0.85 * hh + hh) as i32,
                    (8, 8),
                    &format!(
                        "FPS:{:#0.01} || DT:{:#0.04} || time:{:#0.01}",
                        render_timer.fps, render_timer.delta_time, render_timer.total_time
                    ),
                    Some(&tre::Pixel::new(240, 20, 60, 200)),
                );

                // pixels

                copy_v_frame(pixels.get_frame(), &mut r);

                if pixels.render().is_err() {
                    *control_flow = ControlFlow::Exit;
                    return;
                }
                render_timer.end();
            }
            Event::RedrawEventsCleared => {}
            Event::LoopDestroyed => {}
            _ => {}
        }

        //
    });
}

/// copy's and flips render engines color buffer to pixels buffer
fn copy_v_frame(frame: &mut [u8], re: &mut tre::RenderEngine) {
    for (i, pixel) in frame.chunks_exact_mut(4).enumerate() {
        let (x, y) = tre::index_to_xy(i, BUFF_WIDTH as usize);

        let ny = remap(y as i32, 0, BUFF_HEIGHT as i32, BUFF_HEIGHT as i32, 0) - 1;

        let c = re.buffer.get_pixel(x as u32, ny as u32);

        let rgba = [c.r, c.g, c.b, c.a];

        pixel.copy_from_slice(&rgba);
    }
}
