# rust-fluid-sim
A simple fluid simulation built in rust

This uses a naive implementation where every particle is pushed or pulled towards each other based on which effect is stronger at a given distance.
For this reason, it is not an accurate simulation.

I used a bounding volume hierarchy implementation from [rstar](https://github.com/georust/rstar)
 for performance.

The particles are rendered using [bevy](https://github.com/bevyengine/bevy).
