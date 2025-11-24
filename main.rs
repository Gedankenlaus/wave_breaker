use core::f64;
use log::info;

use iced::mouse;
use iced::widget::canvas;
use iced::widget::canvas::Geometry;
use iced::window;
use iced::{Element, Fill, Rectangle, Renderer, Size, Subscription, Theme};

use num_complex::{Complex64, c64};
use rand::seq::SliceRandom;
use rand::{Rng, random};

use std::io::Cursor;
use std::time::Instant;

use itertools::iproduct;

struct Field2d {
    grid: Vec<Vec<Complex64>>,
    width: usize,
    height: usize,
}

impl Field2d {
    pub fn new(size_x: usize, size_y: usize) -> Field2d {
        Field2d {
            grid: vec![vec![c64(0.0, 0.0); size_x]; size_y],
            width: size_x,
            height: size_y,
        }
    }

    pub fn fill_column(&mut self, column_nr: usize, value: Complex64) {
        for i in 0..self.height {
            self.grid[i][column_nr] = value;
        }
    }

    pub fn double_slit_column(
        &mut self,
        column_nr: usize,
        slit0_0: usize,
        slit0_1: usize,
        slit1_0: usize,
        slit1_1: usize,
    ) {
        for i in 0..slit0_0 {
            self.grid[i][column_nr] = c64(0.0, 0.0);
        }

        for i in slit0_1..slit1_0 {
            self.grid[i][column_nr] = c64(0.0, 0.0);
        }

        for i in slit1_1..self.height {
            self.grid[i][column_nr] = c64(0.0, 0.0);
        }
    }

    pub fn convolve_3x3(&self, stencil: &[[f64; 3]; 3], result_field: &mut Field2d) {
        if result_field.width != self.width || result_field.height != self.height {
            panic!("result_field and this field must have the same dimensions!");
        }
        for y_pos in 0..self.height {
            for x_pos in 0..self.width {
                let mut total_conv_value = c64(0.0, 0.0);

                // Convolve the 3x3 stencil
                for y_stencil in 0..stencil.len() {
                    let mut field_y_index = y_pos + y_stencil;
                    if field_y_index < 1 {
                    } else if field_y_index >= self.height {
                        field_y_index = self.height - 1;
                    } else {
                        field_y_index = field_y_index - 1;
                    }

                    for x_stencil in 0..stencil[y_stencil].len() {
                        let mut field_x_index = x_pos + x_stencil;
                        if field_x_index < 1 {
                            field_x_index = 0;
                        } else if field_x_index >= self.width {
                            field_x_index = self.width - 1;
                        } else {
                            field_x_index = field_x_index - 1;
                        }

                        total_conv_value +=
                            self.grid[field_y_index][field_x_index] * stencil[y_stencil][x_stencil];
                    }
                }
                result_field.grid[y_pos][x_pos] = total_conv_value;
            }
        }
    }
}

fn main() -> iced::Result {
    iced::application(
        "Wave Field - 2D Version",
        WaveField::update,
        WaveField::view,
    )
    .subscription(WaveField::subscription)
    .theme(WaveField::theme)
    .window(window::Settings {
        size: Size {
            width: 1024.,
            height: 512.,
        },
        resizable: false,
        ..window::Settings::default()
    })
    .run()
}

struct WaveField {
    state: State,
}

#[derive(Debug, Clone, Copy)]
enum Message {
    Tick(Instant),
}

impl WaveField {
    fn new() -> Self {
        Self {
            state: State::new(1024, 512),
        }
    }

    fn update(&mut self, message: Message) {
        match message {
            Message::Tick(instant) => {
                self.state.update(instant);
            }
        }
    }

    fn view(&self) -> Element<'_, Message> {
        canvas(&self.state).width(Fill).height(Fill).into()
    }

    fn theme(&self) -> Theme {
        Theme::Moonfly
    }

    fn subscription(&self) -> Subscription<Message> {
        window::frames().map(Message::Tick)
    }
}

impl Default for WaveField {
    fn default() -> Self {
        WaveField::new()
    }
}

struct State {
    start: Instant,
    now: Instant,
    wave_energy_field: Field2d,
    wave_momentum_field: Field2d,
    laplacian_field_buffer: [Field2d; 2],
    li: usize,
    wave_field_cache: canvas::Cache,
    field_image: image::RgbImage,
    field_img_handle: iced::widget::image::Handle,
    total_time: f64,
    field_collpased: bool,
    captured_particles: Vec<[usize; 2]>,
}

impl State {
    const DELTA_T: f64 = 1.0 / 2.0;
    const LAPLACE_STENCIL: [[f64; 3]; 3] = [[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]];
    const LAMBDA: f64 = 3.0;
    const PROP_SPEED: f64 = 10.0;

    pub fn new(field_width: usize, field_height: usize) -> State {
        let now = Instant::now();
        let mut wave_energy_field = Field2d::new(field_width, field_height);
        let mut wave_momentum_field = Field2d::new(field_width, field_height);

        let radius_sq = 10.0 * 10.0;
        let part_area = radius_sq * f64::consts::PI;
        let part_orig_x = 370.0;
        let part_orig_y = 256.0;
        let mut total_prob = 0.0;
        for y in 230..280 {
            for x in 350..390 {
                let transf_x = part_orig_x - (x as f64);
                let transf_y = part_orig_y - (y as f64);

                if (transf_x * transf_x + transf_y * transf_y) <= radius_sq {
                    let wave_number = 2.0 * f64::consts::PI / State::LAMBDA;
                    let arg = Complex64::cis(wave_number * (x as f64)) / part_area;
                    total_prob += 1.0;
                    wave_energy_field.grid[y][x] = arg;
                    wave_momentum_field.grid[y][x] =
                        -1.0 * c64(0.0, wave_number * State::PROP_SPEED) * arg;
                }
            }
        }

        println!("The total probability is: {}", total_prob / part_area);

        // I need an interference pattern at 512 to 522 and i basically want two holes
        // midpoint is 256, a is 60, so at from 196-226 and 286-316
        for i in 512..522 {
            wave_energy_field.double_slit_column(i, 236, 246, 266, 276);
        }

        State {
            start: now,
            now,
            wave_energy_field,
            wave_momentum_field,
            laplacian_field_buffer: [
                Field2d::new(field_width, field_height),
                Field2d::new(field_width, field_height),
            ],
            li: 0,
            wave_field_cache: canvas::Cache::default(),
            field_image: image::RgbImage::new(field_width as u32, field_height as u32),
            field_img_handle: iced::widget::image::Handle::from_rgba(
                1,
                1,
                vec![0u8, 0u8, 0u8, 0u8],
            ),
            total_time: 0.0,
            field_collpased: false,
            captured_particles: vec![],
        }
    }

    fn reset_state(&mut self) {
        self.now = Instant::now();
        let mut wave_energy_field = Field2d::new(1024, 512);
        let mut wave_momentum_field = Field2d::new(1024, 512);

        // for i in 0..512{
        //     //period over 20 pxs
        //     // 2 pi / 20 * i
        //     let arg = Complex64::cis(2.0 * f64::consts::PI / 20.0 * (i as f64));

        //     wave_energy_field.fill_column(i, arg);
        // }

        let radius_sq = 10.0 * 10.0;
        let part_area_inv = 1.0 / (radius_sq * f64::consts::PI);
        let part_orig_x = 370.0;
        let part_orig_y = 256.0;

        for y in 230..280 {
            for x in 350..390 {
                let transf_x = part_orig_x - (x as f64);
                let transf_y = part_orig_y - (y as f64);

                if (transf_x * transf_x + transf_y * transf_y) <= radius_sq {
                    let wave_number = 2.0 * f64::consts::PI / State::LAMBDA;
                    let arg = Complex64::cis(wave_number * (x as f64)) * part_area_inv;
                    wave_energy_field.grid[y][x] = arg;
                    wave_momentum_field.grid[y][x] =
                        -1.0 * c64(0.0, wave_number * State::PROP_SPEED) * arg;
                }
            }
        }

        // I need an interference pattern at 512 to 522 and i basically want two holes
        // midpoint is 256, a is 60, so at from 196-226 and 286-316
        for i in 512..522 {
            wave_energy_field.double_slit_column(i, 236, 246, 266, 276);
        }

        self.start = self.now;
        self.wave_energy_field = wave_energy_field;
        self.wave_momentum_field = wave_momentum_field;
        self.laplacian_field_buffer = [Field2d::new(1024, 512), Field2d::new(1024, 512)];
        self.li = 0;
        self.total_time = 0.0;
        self.field_collpased = false;
    }

    fn collapse_field(&mut self, x_coll: usize, y_coll: usize) {
        let mut wave_energy_field = Field2d::new(1024, 512);
        let mut wave_momentum_field = Field2d::new(1024, 512);

        self.captured_particles.push([x_coll, y_coll]);

        // for i in 0..512{
        //     //period over 20 pxs
        //     // 2 pi / 20 * i
        //     let arg = Complex64::cis(2.0 * f64::consts::PI / 20.0 * (i as f64));

        //     wave_energy_field.fill_column(i, arg);
        // }

        let omega = 2.0 * f64::consts::PI / State::LAMBDA;
        let wave_value = Complex64::cis(omega * (x_coll as f64));

        wave_energy_field.grid[y_coll][x_coll] = wave_value;
        wave_momentum_field.grid[y_coll][x_coll] = -1.0 * c64(0.0, omega) * wave_value;

        // I need an interference pattern at 512 to 522 and i basically want two holes
        // midpoint is 256, a is 60, so at from 196-226 and 286-316
        for i in 512..522 {
            wave_energy_field.double_slit_column(i, 236, 246, 266, 276);
        }

        self.start = self.now;
        self.wave_energy_field = wave_energy_field;
        self.wave_momentum_field = wave_momentum_field;
        self.laplacian_field_buffer = [Field2d::new(1024, 512), Field2d::new(1024, 512)];
        self.li = 0;
        self.field_collpased = true;
    }

    pub fn update(&mut self, now: Instant) {
        self.start = self.start.min(now);
        self.now = now;
        self.wave_field_cache.clear();
        self.total_time += State::DELTA_T;

        if self.total_time >= 600.0 || self.field_collpased {
            self.reset_state();
        }

        // Perform leap frog integration here
        // integrate energy
        for y_index in 0..self.wave_energy_field.height {
            for x_index in 0..self.wave_energy_field.width {
                let current_energy = self.wave_energy_field.grid[y_index][x_index];
                let current_momentum = self.wave_momentum_field.grid[y_index][x_index];
                let current_laplacian = self.laplacian_field_buffer[self.li].grid[y_index][x_index];

                self.wave_energy_field.grid[y_index][x_index] = current_energy
                    + current_momentum * Self::DELTA_T
                    + 0.5 * current_laplacian * Self::DELTA_T * Self::DELTA_T;
            }
        }

        // here, update the first row with a periodic oscilator
        // psi(x,t) = A * exp[i*(2*pi/lambda*x - 2*pi/T*t)]
        // we only evaluate the real part here, in the future it should probably be a complex value :)
        // let psi_value = arg.cos().abs();
        // self.wave_energy_field.fill_column(0, 1.0);

        // reinforce double slit column
        for i in 512..522 {
            self.wave_energy_field
                .double_slit_column(i, 236, 246, 266, 276);
        }

        let li_prev = self.li;

        self.li = (self.li + 1) % self.laplacian_field_buffer.len();

        // compute update second-order deriviative for momentum integration
        self.wave_energy_field.convolve_3x3(
            &Self::LAPLACE_STENCIL,
            &mut self.laplacian_field_buffer[self.li],
        );

        for y_index in 0..self.wave_energy_field.height {
            for x_index in 0..self.wave_energy_field.width {
                let current_momentum = self.wave_momentum_field.grid[y_index][x_index];
                let current_laplacian = self.laplacian_field_buffer[li_prev].grid[y_index][x_index];
                let next_laplacian = self.laplacian_field_buffer[self.li].grid[y_index][x_index];

                self.wave_momentum_field.grid[y_index][x_index] =
                    current_momentum + 0.5 * (current_laplacian + next_laplacian) * Self::DELTA_T;
            }
        }

        let mut rng = rand::rng();

        // detect whether a particle hits in one slit

        let y_range = 0..512;
        let x_range = 712..732;
        let mut total_range: Vec<[usize; 2]> =
            iproduct!(y_range, x_range).map(|(a, b)| [a, b]).collect();
        total_range.shuffle(&mut rng);

        if !self.field_collpased {
            for random_range in total_range.iter() {
                let y_det = random_range[0];
                let x_det = random_range[1];
                let field_value = self.wave_energy_field.grid[y_det][x_det].norm_sqr();
                let random_number = rng.random_range(0.0..1.0);
                if random_number < field_value * 10.0 {
                    self.collapse_field(x_det, y_det);
                    println!("Detected a particle (with probability {field_value})");
                }
            }
        }

        for y_index in 0..self.wave_energy_field.height {
            for x_index in 0..self.wave_energy_field.width {
                let field_value =
                    self.wave_energy_field.grid[y_index][x_index] * 100.0 * f64::consts::PI;
                let field_re = field_value.re;
                let field_im = field_value.im;
                let converted_re_value = 255.0 * field_re.abs();
                let clamped_re_value = converted_re_value.max(0.0).min(255.0);
                let converted_im_value = 255.0 * field_im.abs();
                let clamped_im_value = converted_im_value.max(0.0).min(255.0);

                let mut r_value = 0.0;
                let g_value: f64;
                let mut b_value = 0.0;
                if field_re < 0.0 {
                    r_value = clamped_re_value;
                } else {
                    b_value = clamped_re_value;
                }
                if field_im > 0.0 {
                    g_value = clamped_im_value;
                } else {
                    r_value = r_value + clamped_im_value;
                    g_value = clamped_im_value;
                }

                // r_value /= 2.0;
                // g_value /= 2.0;
                // b_value /= 2.0;

                let rgb_val = image::Rgb([r_value as u8, g_value as u8, b_value as u8]);
                self.field_image
                    .put_pixel(x_index as u32, y_index as u32, rgb_val);
            }
        }

        for random_range in total_range.iter() {
            let x_index = random_range[1] as u32;
            let y_index = random_range[0] as u32;
            let rgb_val = image::Rgb([63u8, 0u8, 127u8]);
            self.field_image.put_pixel(x_index, y_index, rgb_val);
        }

        for collapsed_pix in self.captured_particles.iter() {
            let x_index = collapsed_pix[0] as u32;
            let y_index = collapsed_pix[1] as u32;
            let rgb_val = image::Rgb([255u8, 255u8, 255u8]);
            self.field_image.put_pixel(x_index, y_index, rgb_val);
        }

        for ds_x in 512..522 {
            let rgb_val = image::Rgb([50u8, 50u8, 50u8]);
            for ds_y in 0..236 {
                self.field_image.put_pixel(ds_x, ds_y, rgb_val);
            }

            for ds_y in 245..266 {
                self.field_image.put_pixel(ds_x, ds_y, rgb_val);
            }

            for ds_y in 276..512 {
                self.field_image.put_pixel(ds_x, ds_y, rgb_val);
            }
        }

        let mut cursor = Cursor::new(Vec::new());
        self.field_image
            .write_to(&mut cursor, image::ImageFormat::Png)
            .expect("Failed to encode image data to memory");
        self.field_img_handle = iced::widget::image::Handle::from_bytes(cursor.into_inner());
    }
}

impl<Message> canvas::Program<Message> for State {
    type State = ();

    fn draw(
        &self,
        _state: &Self::State,
        renderer: &Renderer,
        _theme: &Theme,
        bounds: Rectangle,
        _cursor: mouse::Cursor,
    ) -> Vec<Geometry> {
        let frame_buffer = self
            .wave_field_cache
            .draw(renderer, bounds.size(), |frame| {
                //frame.draw_image(bounds, self.field_image);
                frame.draw_image(bounds, &self.field_img_handle);
            });

        vec![frame_buffer]
    }
}
