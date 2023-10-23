//! Simple particle simulation in Rust
use bevy::ecs::bundle;
use bevy::ecs::query::WorldQuery;
use bevy::window::PrimaryWindow;
use bevy::{ecs::query, prelude::*, sprite::MaterialMesh2dBundle};
use nalgebra::Point2;
use rand::prelude::*;
use rayon::prelude::*;
use rstar::{PointDistance, RTree, RTreeObject, AABB};

const PARTICLECOUNT: u32 = 8000; // the number of particles to simulate
const PARTICLEREPULSION: f32 = 0.03; //0.04;
const PARTICLEATTRACTION: f32 = 0.016; // the attraction force between particles aka the y offset of the repulsion function
const PARTICLEMAXREPULSION: f32 = 10.0; // the maximum repulsion force that can be applied to a particle by another particle
const PARTICLETERMINALVELOCITY: f32 = 10.0;
const PARTICLESPHEREOFINFLUENCE: f32 = 30.0; // any particle within this distance will be affected by the repulsion
const PARTICLEREPULSIONDISTANCE: f32 = 0.1; // the width of the graph of the repulsion function
const PARTICLEDRAG: f32 = 0.01; // the amount of force applied to the particle to slow it down

const BOXSIZE: [f32; 2] = [1200., 700.]; // Width, Height
const BOXBOUNCINESS: f32 = 0.9;
const BOXREPULSION: f32 = 10.0;
const BOXMAXREPULSION: f32 = 10.0;
const GRAVITY: f32 = 1.28; // 0.32; // 9.81 / 60.0;

const TIMESTEP: f32 = 1. / 2.; // the amount of time to simulate per frame
const STEPSPERFRAME: u32 = 1;
const DELTATIME: f32 = TIMESTEP / STEPSPERFRAME as f32; // the amount of time to simulate per step

const MOUSEATTRACTFORCE: f32 = 1.0; // the amount of force to attract the particles with

// Create a new type for a particle of water with velocity and position
#[derive(Copy, Clone, Debug, Component)]
struct Particle {
    position: Point2<f32>,
    velocity: Point2<f32>,
    index: u32,
}

// implement RTreeObject for Point2<f32>
impl RTreeObject for Particle {
    type Envelope = AABB<[f32; 2]>;

    fn envelope(&self) -> Self::Envelope {
        AABB::from_corners(
            [
                self.position.x - (2.0_f32.sqrt() * PARTICLESPHEREOFINFLUENCE),
                self.position.y - (2.0_f32.sqrt() * PARTICLESPHEREOFINFLUENCE),
            ],
            [
                self.position.x + (2.0_f32.sqrt() * PARTICLESPHEREOFINFLUENCE),
                self.position.y + (2.0_f32.sqrt() * PARTICLESPHEREOFINFLUENCE),
            ],
        )
    }
}

impl PointDistance for Particle {
    fn distance_2(&self, point: &[f32; 2]) -> f32 {
        let dx = self.position.x - point[0];
        let dy = self.position.y - point[1];
        dx * dx + dy * dy // returning squared distance is required
    }

    fn contains_point(&self, point: &[f32; 2]) -> bool {
        self.distance_2(point) <= PARTICLESPHEREOFINFLUENCE * PARTICLESPHEREOFINFLUENCE
    }
}

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
                position: Point2::new(0.0, 0.0),
                velocity: Point2::new(0.0, 0.0),
                index: i,
            });
    }
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_systems(PreStartup, setup_circles)
        .add_systems(Startup, initialize_positions)
        .add_systems(Update, update)
        .run();
}

fn initialize_positions(mut particleQuery: Query<&mut Particle, With<Particle>>) {
    let mut index = 0;
    for x in 0..((PARTICLECOUNT as f32).sqrt() as u16) {
        for y in 0..((PARTICLECOUNT as f32).sqrt() as u16) {
            let mut particle = particleQuery.iter_mut().nth(index).unwrap(); // extremely inefficient but only run once
            particle.position.x = (((x * 2) as f32) - 500.).into();
            particle.position.y = (((y * 2) as f32) - 250.).into();
            index += 1;
        }
    }
}

fn update(
    mut transfrom_query: Query<&mut Transform, With<Particle>>,
    mut particle_query: Query<&mut Particle>,
    window_query: Query<&Window, With<PrimaryWindow>>,
    buttons: Res<Input<MouseButton>>,
)
// Called every frame
{
    let mut particles = particle_query.iter_mut().collect::<Vec<_>>();

    let transforms = transfrom_query.iter_mut().collect::<Vec<_>>();

    for _ in 0..STEPSPERFRAME {
        let mut position = Vec2::new(0., 0.);
        if buttons.pressed(MouseButton::Left) || buttons.pressed(MouseButton::Right) {
            // Left Button is being held down
            if let Some(cursor_position) = window_query.single().cursor_position() {
                let mut cursor_position_2 = cursor_position + Vec2::new(-600., -350.);
                cursor_position_2.y = -cursor_position_2.y;
                position = cursor_position_2;
            }
        }

        // make a copy of the particles to read from in Vec<Particle>
        let read_only_particles: Vec<Particle> = particles
            .iter()
            .map(|particle| Particle {
                position: particle.position,
                velocity: particle.velocity,
                index: particle.index,
            })
            .collect();

        // format the particles into an RTree
        let particle_tree = RTree::bulk_load(read_only_particles);

        particles.par_iter_mut().for_each(|particle| {
            simulation_step(
                particle,
                &particle_tree,
                &position,
                buttons.pressed(MouseButton::Left),
            );
        });
    }

    set_circle_positions(particles, transforms); // only time transforms is updated meaning steps per frame doesn't work.
}

fn set_circle_positions(
    mut particles: Vec<Mut<'_, Particle>>,
    mut circles: Vec<Mut<'_, Transform>>,
) {
    particles.iter_mut().enumerate().for_each(|(i, particle)| {
        circles[i as usize].translation.x = particle.position.x;
        circles[i as usize].translation.y = particle.position.y;
    });
}

fn simulation_step(
    particle: &mut Particle,
    particle_tree: &RTree<Particle>,
    cursor_position: &Vec2,
    is_left_button_pressed: bool,
) {
    drag_step(particle);
    gravity_step(particle);
    cursor_attraction_step(particle, cursor_position, is_left_button_pressed);
    particle_repulsion_step(particle, particle_tree);
    box_repulsion_step(particle);
    bounce_step(particle);
    cap_velocity(particle);
    apply_velocity_step(particle);
}

fn drag_step(particle: &mut Particle) {
    particle.velocity.x *= 1.0 - (PARTICLEDRAG * DELTATIME);
    particle.velocity.y *= 1.0 - (PARTICLEDRAG * DELTATIME);
}

fn gravity_step(particle: &mut Particle) {
    particle.velocity.y -= GRAVITY * DELTATIME;
}

fn bounce_step(particle: &mut Particle) {
    if particle.position.x.abs() > BOXSIZE[0] / 2. {
        particle.velocity.x = (-particle.velocity.x) * BOXBOUNCINESS;
    }
    if particle.position.y.abs() > BOXSIZE[1] / 2. {
        particle.velocity.y = (-particle.velocity.y) * BOXBOUNCINESS;
    }
}

fn box_repulsion_step(particle: &mut Particle) {
    // check if the particle is outside the box
    if particle.position.x.abs() > BOXSIZE[0] / 2. || particle.position.y.abs() > BOXSIZE[1] / 2. {
        return;
    }

    let nearest_rectangle_point =
        nearest_point_on_rectangle([particle.position.x, particle.position.y], BOXSIZE);

    let dx = particle.position.x - nearest_rectangle_point[0];
    let dy = particle.position.y - nearest_rectangle_point[1];

    let distance = (dx * dx + dy * dy).sqrt();

    particle.velocity.x +=
        (dx / distance.powi(2)).min(BOXMAXREPULSION) * BOXREPULSION as f32 * DELTATIME;
    particle.velocity.y +=
        (dy / distance.powi(2)).min(BOXMAXREPULSION) * BOXREPULSION as f32 * DELTATIME;
}

#[inline(never)]
fn particle_repulsion_step(particle: &mut Particle, particle_tree: &RTree<Particle>) {
    let mut total_force = [0.0, 0.0];
    let repulsion_distance_squared = PARTICLEREPULSIONDISTANCE.powi(2);
    let particles = particle_tree.locate_all_at_point(&[particle.position.x, particle.position.y]);
    for other_particle in particles {
        // make sure the particle is not itself
        if particle.index != other_particle.index {
            let dx = particle.position.x - other_particle.position.x;
            let dy = particle.position.y - other_particle.position.y;

            let distance_squared = dx.powi(2) + dy.powi(2);

            if distance_squared < 0.05 {
                total_force[0] += rand::thread_rng().gen::<f32>() * PARTICLEMAXREPULSION;
                total_force[1] += rand::thread_rng().gen::<f32>() * PARTICLEMAXREPULSION;
            } else {
                let force = ((1.0 / (distance_squared * repulsion_distance_squared))
                    - PARTICLEATTRACTION)
                    .min(PARTICLEMAXREPULSION);

                total_force[0] += dx * force;
                total_force[1] += dy * force;
            }
        }
    }
    particle.velocity.x += total_force[0] * PARTICLEREPULSION * DELTATIME;
    particle.velocity.y += total_force[1] * PARTICLEREPULSION * DELTATIME;
}

fn cursor_attraction_step(
    particle: &mut Particle,
    cursor_position: &Vec2,
    is_left_button_pressed: bool,
) {
    if is_left_button_pressed && *cursor_position != Vec2::new(0., 0.) {
        // Left Button is being held down
        let dx = particle.position.x - cursor_position.x as f32;
        let dy = particle.position.y - cursor_position.y as f32;

        let distance_squared = dx.powi(2) + dy.powi(2);

        let force = (distance_squared * MOUSEATTRACTFORCE.powi(2)).min(10000.0);

        particle.velocity.x -= dx * force * DELTATIME;
        particle.velocity.y -= dy * force * DELTATIME;
    } else if *cursor_position != Vec2::new(0., 0.) {
        // Right Button is being held down
        let dx = particle.position.x - cursor_position.x as f32;
        let dy = particle.position.y - cursor_position.y as f32;

        let distance_squared = dx.powi(2) + dy.powi(2);

        let force = (distance_squared * MOUSEATTRACTFORCE.powi(2)).min(10000.0);

        particle.velocity.x += dx * force * DELTATIME;
        particle.velocity.y += dy * force * DELTATIME;
    }
}

fn apply_velocity_step(particle: &mut Particle) {
    particle.position.x += particle.velocity.x * DELTATIME;
    particle.position.y += particle.velocity.y * DELTATIME;
}

fn cap_velocity(particle: &mut Particle) {
    // get the magnitude of the velocity
    let velocity_magnitude = (particle.velocity.x.powi(2) + particle.velocity.y.powi(2)).sqrt();
    if velocity_magnitude > PARTICLETERMINALVELOCITY {
        particle.velocity.x *= PARTICLETERMINALVELOCITY / velocity_magnitude;
        particle.velocity.y *= PARTICLETERMINALVELOCITY / velocity_magnitude;
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
