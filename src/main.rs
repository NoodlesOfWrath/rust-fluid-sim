//! Shows how to render simple primitive shapes with a single color.

const PARTICLECOUNT: u32 = 1000;

use bevy::{prelude::*, sprite::MaterialMesh2dBundle};
use std::num;

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
) {
    commands.spawn(Camera2dBundle::default());

    // Make an array to store all circles

    for _ in 0..PARTICLECOUNT {
        // Circle
        let mut circle = commands.spawn(MaterialMesh2dBundle {
            mesh: meshes.add(shape::Circle::new(5.).into()).into(),
            material: materials.add(ColorMaterial::from(Color::BLUE)),
            transform: Transform::from_translation(Vec3::new(0., 0., 0.)),
            ..default()
        });
    }
}

// Create a new type for a particle of water with velocity and position
#[derive(Copy, Clone, Debug)]
struct Particle {
    position: (f32, f32),
    velocity: (f32, f32),
}

const BOXDIMS: [i32; 4] = [-10, 10, -10, 10]; // Left, Right, Bottom, Top

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_systems(Startup, setup)
        .run();

    let mut particles: [Particle; PARTICLECOUNT as usize] = [Particle {
        position: (0.0, 0.0),
        velocity: (0.0, 0.0),
    }; 1000];

    for i in 0..PARTICLECOUNT {
        particles[i as usize] = Particle {
            position: (0.0, 0.0),
            velocity: (2.0, 3.0),
        };
    }
    for mut particle in particles {
        simulation_step(&mut particle);
        velocity_step(&mut particle);
    }
}

fn initialize_positions(particles: &mut [Particle]) {
    let mut index = 0;
    for x in 0..((PARTICLECOUNT as f32).sqrt() as u16) {
        for y in 0..((PARTICLECOUNT as f32).sqrt() as u16) {
            particles[index].position.0 = (((x * 20) as f32) - 500.).into();
            particles[index].position.0 = (((y * 20) as f32) - 250.).into();
            index += 1;
        }
    }
}

fn setCirclePosition() {}

fn simulation_step(particle: &mut Particle) {
    if BOXDIMS[0] as f32 > particle.position.0 || particle.position.0 > BOXDIMS[1] as f32 {
        particle.velocity.0 = -particle.velocity.0;
    }
    if (BOXDIMS[2] as f32 > particle.position.1 || particle.position.1 > BOXDIMS[3] as f32) {
        particle.velocity.1 = -particle.velocity.1;
    }
    velocity_step(particle);
}

fn velocity_step(particle: &mut Particle) {
    particle.position.0 += particle.velocity.0;
    particle.position.1 += particle.velocity.1;
}
