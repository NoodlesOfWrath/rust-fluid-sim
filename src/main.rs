//! Simple particle simulation in Rust

use bevy::ecs::bundle;
use bevy::ecs::query::WorldQuery;
use bevy::{ecs::query, prelude::*, sprite::MaterialMesh2dBundle};
use rand::prelude::*;
use rayon::prelude::*;
use std::num;
use std::time::Instant;

const PARTICLECOUNT: u32 = 8000; // the number of particles to simulate
const PARTICLEREPULSION: f32 = 0.15;
const PARTICLEATTRACTION: f32 = 0.016; // the attraction force between particles aka the y offset of the repulsion function
const PARTICLEMAXREPULSION: f32 = 10.0; // the maximum repulsion force that can be applied to a particle by another particle
const PARTICLETERMINALVELOCITY: f32 = 10.0;
const PARTICLESPHEREOFINFLUENCE: f32 = 30.0; // any particle within this distance will be affected by the repulsion
const PARTICLEREPULSIONDISTANCE: f32 = 0.1; // the width of the graph of the repulsion function
const PARTICLEDRAG: f32 = 0.01; // the amount of force applied to the particle to slow it down

const BOXSIZE: [f32; 2] = [1200., 700.]; // Width, Height
const BOXBOUNCINESS: f32 = 0.9;
const BOXREPULSION: f32 = 1.0;
const BOXMAXREPULSION: f32 = 2.0;
const GRAVITY: f32 = 1.28; // 0.32; // 9.81 / 60.0;

const BOXSEGMENTS: usize = 30; // the number of boxes in each dimension

// Create a new type for a particle of water with velocity and position
#[derive(Copy, Clone, Debug, Component)]
struct Particle {
    position: (f32, f32),
    velocity: (f32, f32),
}

#[derive(Clone, Component)]
struct ParticleBoxes([[Vec<usize>; BOXSEGMENTS]; BOXSEGMENTS]); // Stores the indices of the particles in the boxes

fn setup_circles(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
) {
    commands.spawn(Camera2dBundle::default());

    // Quad
    commands.spawn(MaterialMesh2dBundle {
        mesh: meshes
            .add(shape::Quad::new(Vec2::new(BOXSIZE[0], BOXSIZE[1])).into())
            .into(),
        material: materials.add(ColorMaterial::from(Color::LIME_GREEN)),
        transform: Transform::from_translation(Vec3::new(0., 0., -1.)),
        ..default()
    });

    for i in 1..PARTICLECOUNT {
        // Circle
        commands
            .spawn(MaterialMesh2dBundle {
                mesh: meshes.add(shape::Circle::new(2.).into()).into(),
                material: materials.add(ColorMaterial::from(Color::BLUE)),
                transform: Transform::from_translation(Vec3::new(0., 0., 0.)),
                ..default()
            })
            .insert(Particle {
                position: (0., 0.),
                velocity: (0., 0.),
            });
    }

    // make an object to store the particle boxes
    commands.spawn(ParticleBoxes(std::array::from_fn(|_| {
        std::array::from_fn(|_| vec![])
    })));
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_systems(PreStartup, setup_circles)
        .add_systems(Startup, initialize_positions)
        .add_systems(FixedUpdate, update)
        .run();
}

fn initialize_positions(mut particleQuery: Query<&mut Particle, With<Particle>>) {
    let mut index = 0;
    for x in 0..((PARTICLECOUNT as f32).sqrt() as u16) {
        for y in 0..((PARTICLECOUNT as f32).sqrt() as u16) {
            let mut particle = particleQuery.iter_mut().nth(index).unwrap(); // extremely inefficient but only run once
            particle.position.0 = (((x * 5) as f32) - 500.).into();
            particle.position.1 = (((y * 5) as f32) - 250.).into();
            index += 1;
        }
    }
}

fn update(
    mut transfrom_query: Query<&mut Transform, With<Particle>>,
    mut particle_query: Query<&mut Particle>,
    mut particle_boxes: Query<&mut ParticleBoxes>,
)
// Called every frame
{
    let mut particles = particle_query.iter_mut().collect::<Vec<_>>();
    let mut transforms = transfrom_query.iter_mut().collect::<Vec<_>>();
    let mut particle_boxes = particle_boxes
        .iter_mut()
        .next()
        .map(|x| x.0.clone())
        .unwrap();

    sort_into_boxes(&transforms, &mut particle_boxes);

    particles
        .par_iter_mut()
        .enumerate()
        .for_each(|(i, particle)| {
            simulation_step(particle, &transforms, &particle_boxes, i);
        });

    set_circle_positions(particles, transforms);
}

fn set_circle_positions(
    mut particles: Vec<Mut<'_, Particle>>,
    mut circles: Vec<Mut<'_, Transform>>,
) {
    particles.iter_mut().enumerate().for_each(|(i, particle)| {
        circles[i as usize].translation.x = particle.position.0;
        circles[i as usize].translation.y = particle.position.1;
    });
}

fn simulation_step(
    particle: &mut Particle,
    other_particles: &Vec<Mut<'_, Transform>>,
    particle_boxes: &[[Vec<usize>; BOXSEGMENTS]; BOXSEGMENTS],
    index: usize,
) {
    drag_step(particle);
    gravity_step(particle);
    particle_repulsion_step(particle, other_particles, index, particle_boxes);
    box_repulsion_step(particle);
    bounce_step(particle);
    cap_velocity(particle);
    apply_velocity_step(particle);
}

fn drag_step(particle: &mut Particle) {
    particle.velocity.0 *= 1.0 - PARTICLEDRAG;
    particle.velocity.1 *= 1.0 - PARTICLEDRAG;
}

fn gravity_step(particle: &mut Particle) {
    particle.velocity.1 -= GRAVITY;
}

fn bounce_step(particle: &mut Particle) {
    if particle.position.0.abs() > BOXSIZE[0] / 2. {
        particle.velocity.0 = (-particle.velocity.0) * BOXBOUNCINESS;
    }
    if particle.position.1.abs() > BOXSIZE[1] / 2. {
        particle.velocity.1 = (-particle.velocity.1) * BOXBOUNCINESS;
    }
}

fn box_repulsion_step(particle: &mut Particle) {
    // check if the particle is outside the box
    if particle.position.0.abs() > BOXSIZE[0] / 2. || particle.position.1.abs() > BOXSIZE[1] / 2. {
        return;
    }

    let nearest_rectangle_point =
        nearest_point_on_rectangle([particle.position.0, particle.position.1], BOXSIZE);

    let dx = particle.position.0 - nearest_rectangle_point[0];
    let dy = particle.position.1 - nearest_rectangle_point[1];

    let distance = (dx * dx + dy * dy).sqrt();

    particle.velocity.0 += (dx / distance.powi(2)).min(BOXMAXREPULSION) * BOXREPULSION as f32;
    particle.velocity.1 += (dy / distance.powi(2)).min(BOXMAXREPULSION) * BOXREPULSION as f32;
}

fn particle_repulsion_step(
    particle: &mut Particle,
    other_particles: &Vec<Mut<'_, Transform>>,
    index: usize,
    particle_boxes: &[[Vec<usize>; BOXSEGMENTS]; BOXSEGMENTS],
) {
    // find the 5 boxes that are around the particle
    let [box_x, box_y] = get_box_indices(particle.position.0, particle.position.1);
    let boxes_to_check = [
        [box_x, box_y],
        [box_x + 1, box_y],
        [box_x - 1, box_y],
        [box_x, box_y + 1],
        [box_x, box_y - 1],
    ];

    boxes_to_check.iter().for_each(|[box_x, box_y]| {
        if *box_x > BOXSEGMENTS as i32 - 1
            || *box_x < 0
            || *box_y > BOXSEGMENTS as i32 - 1
            || *box_y < 0
        {
            // if the index is outside the box don't do anything
        } else {
            particle_boxes[*box_x as usize][*box_y as usize]
                .iter()
                .for_each(|other_particle_index| {
                    if *other_particle_index != index {
                        let other_particle = other_particles[*other_particle_index].clone(); // could probably be optimized
                        let dx = particle.position.0 - other_particle.translation.x;
                        let dy = particle.position.1 - other_particle.translation.y;

                        let distance = (dx * dx + dy * dy).sqrt();

                        // Check if the particle is not itself
                        if distance < PARTICLESPHEREOFINFLUENCE {
                            if distance < 0.01 {
                                particle.velocity.0 +=
                                    rand::thread_rng().gen::<f32>() * PARTICLEMAXREPULSION as f32;
                                particle.velocity.1 +=
                                    rand::thread_rng().gen::<f32>() * PARTICLEMAXREPULSION as f32;
                            } else {
                                let force = ((1.0
                                    / (distance * PARTICLEREPULSIONDISTANCE).powi(2))
                                    - PARTICLEATTRACTION)
                                    .min(PARTICLEMAXREPULSION)
                                    * PARTICLEREPULSION as f32;
                                particle.velocity.0 += dx * force;
                                particle.velocity.1 += dy * force;
                            }
                        }
                    }
                });
        }
    });
}

fn apply_velocity_step(particle: &mut Particle) {
    particle.position.0 += particle.velocity.0;
    particle.position.1 += particle.velocity.1;
}

fn cap_velocity(particle: &mut Particle) {
    // get the magnitude of the velocity
    let velocity_magnitude = (particle.velocity.0.powi(2) + particle.velocity.1.powi(2)).sqrt();
    if velocity_magnitude > PARTICLETERMINALVELOCITY {
        particle.velocity.0 *= PARTICLETERMINALVELOCITY / velocity_magnitude;
        particle.velocity.1 *= PARTICLETERMINALVELOCITY / velocity_magnitude;
    }
}

fn nearest_point_on_rectangle(point: [f32; 2], rectangle_size: [f32; 2]) -> [f32; 2] {
    let half_width: f32 = rectangle_size[0] / 2 as f32;
    let half_height: f32 = rectangle_size[1] / 2 as f32;

    // Calculate distances to each edge
    let left_distance = (point[0] + half_width).abs();
    let right_distance = (point[0] - half_width).abs();
    let top_distance = (point[1] + half_height).abs();
    let bottom_distance = (point[1] - half_height).abs();

    // Find the minimum of the four distances
    let min_distance = left_distance
        .min(right_distance)
        .min(top_distance)
        .min(bottom_distance);

    // Return the nearest point and the distance
    if min_distance == left_distance {
        return [-half_width, point[1]];
    } else if min_distance == right_distance {
        return [half_width, point[1]];
    } else if min_distance == top_distance {
        return [point[0], -half_height];
    } else {
        return [point[0], half_height];
    }
}

fn sort_into_boxes(
    particles: &Vec<Mut<'_, Transform>>,
    particle_boxes: &mut [[Vec<usize>; BOXSEGMENTS]; BOXSEGMENTS],
) {
    clear_boxes(particle_boxes);

    // had to use a for loop because the rayon library doesn't support iterators with indices
    particles
        .iter()
        .enumerate()
        .for_each(|(particle_index, particle)| {
            let x = particle.translation.x;
            let y = particle.translation.y;

            let [x_index, y_index] = get_box_indices(x, y);

            particle_boxes[x_index as usize][y_index as usize].push(particle_index);
        });
}

fn get_box_indices(x: f32, y: f32) -> [i32; 2] {
    let x_index = ((x as i32 / (BOXSIZE[0] / (BOXSEGMENTS / 2) as f32) as i32)
        + BOXSEGMENTS as i32 / 2)
        .min(BOXSEGMENTS as i32 - 1)
        .max(0);
    let y_index = ((y as i32 / (BOXSIZE[1] / (BOXSEGMENTS / 2) as f32) as i32)
        + BOXSEGMENTS as i32 / 2)
        .min(BOXSEGMENTS as i32 - 1)
        .max(0);
    return [x_index, y_index];
}

fn clear_boxes(particle_boxes: &mut [[Vec<usize>; BOXSEGMENTS]; BOXSEGMENTS]) {
    particle_boxes.iter_mut().for_each(|x| {
        x.par_iter_mut().for_each(|y| {
            y.clear();
        });
    });
}
