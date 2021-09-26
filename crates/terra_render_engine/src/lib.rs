mod camera;
mod mesh;
mod pixel;
mod ray_tracing;
mod render_engine;
mod retro_bit_map;
mod utility;

// what to show, dumps content of mod in lib
pub use camera::*;
pub use mesh::*;
pub use pixel::*;
pub use ray_tracing::*;
pub use render_engine::*;
pub use retro_bit_map::*;
pub use utility::*;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
