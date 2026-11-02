#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//extensions: Colorful Comments, indent-rainbow,  Party Owl '84, vscode-pets

// ! IMPORTANT: All team members must use SI units consistently!
// !
// ! Position / Distance:   meters (m)
// ! Velocity:              meters per second (m/s)
// ! Acceleration:          meters per second squared (m/s^2)
// ! Time / Time step:      seconds (s)
// ! Mass:                  kilograms (kg)
// ! Electric charge:       coulombs (C)
// ! Force:                 newtons (N)
// ! Particles              Array that stores the data of the particles
// ! All physical calculations must use these units to ensure consistency between the different simulation components.
// ! Note: Image dimensions are measured in pixels (px), and
// ! atomic number (Z) is dimensionless.


//-------------------------
//* Definition of structures
//-------------------------

//can be used to store position(x, y), velocity(vx, vy), force(Fx, Fy), accelaration(ax,ay)
typedef struct
{
    double x;
    double y;
} Vector2;
//position, atomic number (positive charge, defaults to 79 if missing)
typedef struct {
    double x;
    double y;
    int Z;
} Nucleus;
//color definition for unique colors
typedef struct {
    unsigned char r;
    unsigned char g;
    unsigned char b;
} RGB;
//stores data for each alpha particle
typedef struct {
    Vector2 position;
    Vector2 velocity;
    RGB COLOR;
} AlphaParticle;

//-------------------------
//* Definition of constants
//-------------------------

#define COULOMB_CONSTANT 8.9875517923e9
#define ELEMENTARY_CHARGE 1.602176634e-19
#define ALPHA_CHARGE (2.0 * ELEMENTARY_CHARGE)
#define ALPHA_MASS 6.644657e-27

//-------------------------
//* Functions of phase 1
//-------------------------

//calculates the coulomb force btw an alpha particle and a nucleus. The position is the position of an Alpha particle, but only its position should be passed in.
Vector2 calculate_coulomb_force(Vector2 position, Nucleus nucleus){
    double Fx, Fy, dx, dy, d;
    Vector2 v;

    dx = position.x - nucleus.x;
    dy = position.y - nucleus.y;
    d = sqrt(dx * dx + dy * dy);
    //the columb force multiplied by the direction components, sin(dy/d) or the cos(dx/d), but bc we dont know the angle we use the dx formula
    Fx = COULOMB_CONSTANT * (ALPHA_CHARGE * nucleus.Z * ELEMENTARY_CHARGE) * dx / (d * d * d);
    Fy = COULOMB_CONSTANT * (ALPHA_CHARGE * nucleus.Z * ELEMENTARY_CHARGE) * dy / (d * d * d);

    v.x = Fx;
    v.y = Fy;

    return v;
}
//todo: function cant handle very small distances btw particles becouse very small nonzero d can lead to division by 0

Vector2 calculate_total_force(Vector2 position, const Nucleus nuclei[], int nucleus_count){

}

Vector2 calculate_acceleration(Vector2 position, const Nucleus nuclei[], int nucleus_count) {

}

AlphaParticle *initialize_particles(int particle_count,Vector2 *position,Vector2 *velocity) {
    AlphaParticle *Particles;
    int i;
    if (particle_count <= 0 ||position == NULL ||velocity == NULL) {
        return NULL;
        }
    Particles = (AlphaParticle *)malloc(sizeof(AlphaParticle) * particle_count);
    if (Particles == NULL) {
        return NULL;
    }
    for (i = 0; i < particle_count; i++) {
        Particles[i].position.x = position[i].x;
        Particles[i].velocity.x = velocity[i].x; //Zsolt fell free to change to name of this but pls make sure
        Particles[i].velocity.y=0.0;                //to pay attention that here I assume initial_velocity in a way such that it has only x component
        Particles[i].position.y = (i + 0.5) * physical_height / particle_count; //We must decide on the height of the playing field (here I used it as double)
        Particles[i].COLOR.r =rand() % 256;
        Particles[i].COLOR.g =rand() % 256;
        Particles[i].COLOR.b =rand() % 256;
    }
    return Particles;
}
AlphaParticle *simulate_particles(Nucleus Nuclei, AlphaParticle *Particles, int particle_count) {

}
//-------------------------
//* Functions of phase 2
//-------------------------


//-------------------------
//* Functions of phase 3
//-------------------------


//-------------------------
//* Main
//-------------------------
int main() {
    int nucleus_count = 0;// will be scanned from file, 0 is only placeholder
    Nucleus nuclei[nucleus_count];
    AlphaParticle *Particles; // Array that will store Particles data
    // Phase 2: User Input & File Handling
    read_parameters();       // Zsolt
    load_nuclei();           // Ábel

    // Phase 1: Physics & Simulation
    initialize_particles();  // Pityesz
    simulate_particles();    // Pityesz
    //&Inside simulate_particles(), Pityesz will call Ábel's rk4_step(), which will use Zsolt's three physics functions:

    // Phase 3: Visualization & Output
    create_image();          // Zsolt
    draw_nuclei();           // Ábel
    draw_trajectories();     // Ábel
    save_ppm();              // Pityesz

    // Clean up
    free_memory();           // Shared responsibility

    return  0;
}