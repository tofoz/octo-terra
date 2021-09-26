use std::time::Instant;

pub struct FrameTimer {
    pub total_time: f64,
    pub delta_time: f64,
    pub fps: f64,
    smoothing: f64,
    fps_timer: Instant,
}

impl FrameTimer {
    pub fn new(smoothing: f64) -> Self {
        Self {
            total_time: 0.0,
            delta_time: 0.00001,
            fps: 0.0,
            smoothing,
            fps_timer: Instant::now(),
        }
    }
    pub fn start(&mut self) {
        self.fps_timer = Instant::now();
    }
    pub fn end(&mut self) {
        if self.total_time >= f64::MAX - 1.0 {
            self.total_time = 0.0;
        }
        self.total_time += self.delta_time;
        self.delta_time = self.fps_timer.elapsed().as_secs_f64();

        self.fps = (1. / self.delta_time * self.smoothing) + (self.fps * (1.0 - self.smoothing));
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}

// windowing
trait State {
    fn resize(&mut self) {}
    fn input() {}
    fn update() {}
    fn render() {}
}
