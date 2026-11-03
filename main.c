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
// !
// ! All physical calculations must use these units to ensure consistency between the different simulation components.
// ! Note: Image dimensions are measured in pixels (px), and
// ! atomic number (Z) is dimensionless.

// ? Parameter	Value
// ? Image width	800 pixels
// ? Image height	800 pixels
// ? Physical width	1e-9 m
// ? Physical height	1e-9 m
// ? Starting X	0 m
// ? Initial velocity	1e7 m/s
// ? Time step (dt)	1e-19 s
// ? Maximum steps	1000


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
#define PI 3.14159265358979
#define EPSILON_0 8.8541878188e-12
#define COULOMB_CONSTANT (1.0 / (4.0 * PI * EPSILON_0))

#define ELEMENTARY_CHARGE 1.602176634e-19
#define ALPHA_CHARGE (2.0 * ELEMENTARY_CHARGE)
#define ALPHA_MASS 6.644657e-27

//-------------------------
//* Functions of phase 1
//-------------------------

//calculates the coulomb force btw one alpha particle and one nucleus.
//The position is the position of an Alpha particle.
//Returns the total Coulomb force as a Vector2 (Fx, Fy), measured in newtons.
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

Vector2 calculate_total_force(Vector2 position, const Nucleus nuclei[], int nucleus_count) {
    //F(r) = Σi Fi(r)
    Vector2 Sum_of_forces = {0,0};
    Vector2 f;

    for (int i = 0; i < nucleus_count; i++) {
        f = calculate_coulomb_force(position, nuclei[i]);
        Sum_of_forces.x += f.x;
        Sum_of_forces.y += f.y;
    }
    return Sum_of_forces;
}


/* Vector2 calculate_acceleration(Vector2 position, const Nucleus nuclei[], int nucleus_count) {
    Vector2 accelaration_of_particle = {0,0};
    Vector2 F = calculate_total_force(position, nuclei, nucleus_count);
    accelaration_of_particle.x = F.x / ALPHA_MASS;
    accelaration_of_particle.y = F.y / ALPHA_MASS;
    return accelaration_of_particle;
} */


void numerical_calculations(AlphaParticle *alpha, int max_steps, const Nucleus nuclei[], int nucleus_count, double dt){

    int t = 1;
    while (t <= max_steps){

        // first step
        Vector2 k1_r = (*alpha).velocity;   //k1_r = v

        Vector2 k1_v;
        Vector2 F1 = calculate_total_force((*alpha).position, nuclei, nucleus_count);
        k1_v.x = F1.x/ALPHA_MASS;  //k1_v = F(r)/m
        k1_v.y = F1.y/ALPHA_MASS;  //k1_v = F(r)/m


        // second step
        Vector2 k2_r;
        k2_r.x = (*alpha).velocity.x + 0.5 * dt * k1_v.x;  //k2_r = v + 0.5 * dt * k1_v
        k2_r.y = (*alpha).velocity.y + 0.5 * dt * k1_v.y;  //k2_r = v + 0.5 * dt * k1_v

        Vector2 k2_v;
        Vector2 pos2;
        pos2.x = (*alpha).position.x + 0.5*dt*k1_r.x;          //X r + 0.5 * dt * k1_r
        pos2.y = (*alpha).position.y + 0.5*dt*k1_r.y;          //Y r + 0.5 * dt * k1_r
        Vector2 F2 = calculate_total_force((*alpha).position, nuclei, nucleus_count);        //  F(r + 0.5 * dt * k1_r)
        k2_v.x =  F2.x / ALPHA_MASS;                              //X k2_v = F(r + 0.5 * dt * k1_r) / m
        k2_v.y =  F2.y / ALPHA_MASS;                              //Y k2_v = F(r + 0.5 * dt * k1_r) / m


        // third step
        Vector2 k3_r;
        k3_r.x = (*alpha).velocity.x + 0.5*dt*k2_v.x; //X k3_r = v + 0.5 * dt * k2_v
        k3_r.y = (*alpha).velocity.y + 0.5*dt*k2_v.y; //Y k3_r = v + 0.5 * dt * k2_v

        Vector2 k3_v;
        Vector2 pos3;
        pos3.x = (*alpha).position.x + 0.5*dt*k2_r.x;          //X r + 0.5 * dt * k2_r
        pos3.y = (*alpha).position.y + 0.5*dt*k2_r.y;          //Y r + 0.5 * dt * k2_r
        Vector2 F3 = calculate_total_force((*alpha).position, nuclei, nucleus_count);        //  F(r + 0.5 * dt * k2_r)
        k3_v.x =  F3.x / ALPHA_MASS;                              //X k3_v = F(r + 0.5 * dt * k2_r) / m
        k3_v.y =  F3.y / ALPHA_MASS;                              //Y k3_v = F(r + 0.5 * dt * k2_r) / m


        // fourth step
        Vector2 k4_r;
        k4_r.x = (*alpha).velocity.x + dt*k3_v.x; //X k4_r = v + dt * k3_v
        k4_r.y = (*alpha).velocity.y + dt*k3_v.y; //Y k4_r = v + dt * k3_v

        Vector2 k4_v;
        Vector2 pos4;
        pos4.x = (*alpha).position.x + dt*k3_r.x;          //X r + dt * k3_r
        pos4.y = (*alpha).position.y + dt*k3_r.y;          //Y r + dt * k3_r
        Vector2 F4 = calculate_total_force((*alpha).position, nuclei, nucleus_count);        //  F(r + dt * k3_r)
        k4_v.x =  F4.x / ALPHA_MASS;                              //X k4_v = F(r + dt * k3_r) / m
        k4_v.y =  F4.y / ALPHA_MASS;                              //Y k4_v = F(r + dt * k3_r) / m



        //Final update
        Vector2 r_new;
        r_new.x = (*alpha).position.x + dt * (k1_r.x + 2*k2_r.x + 2*k3_r.x + k4_r.x) / 6;  //X r_new = r + dt * (k1_r + 2*k2_r + 2*k3_r + k4_r) / 6
        r_new.y = (*alpha).position.y + dt * (k1_r.y + 2*k2_r.y + 2*k3_r.y + k4_r.y) / 6;  //Y r_new = r + dt * (k1_r + 2*k2_r + 2*k3_r + k4_r) / 6

        Vector2 v_new;
        v_new.x = (*alpha).velocity.x + dt * (k1_v.x + 2*k2_v.x + 2*k3_v.x + k4_v.x) / 6;  //X v_new = v + dt * (k1_v + 2*k2_v + 2*k3_v + k4_v) / 6
        v_new.y = (*alpha).velocity.y + dt * (k1_v.y + 2*k2_v.y + 2*k3_v.y + k4_v.y) / 6;  //Y v_new = v + dt * (k1_v + 2*k2_v + 2*k3_v + k4_v) / 6

        //SAVE IT INTO A FILE (PPM?)

        draw_trajectories((*alpha).position, r_new);

        //increase step
        (*alpha).position = r_new;
        (*alpha).velocity = v_new;
        t++;
    }

    return;
}


AlphaParticle *initialize_particles(int particle_count, Vector2 *position, Vector2 *velocity, int WINDOW_HEIGHT) {
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
        Particles[i].velocity.y = 0.0;                //to pay attention that here I assume initial_velocity in a way such that it has only x component
        Particles[i].position.y = (i + 0.5) * WINDOW_HEIGHT / particle_count; //We must decide on the height of the playing field (here I used it as double)
        Particles[i].COLOR.r =(i + 10) % 256;
        Particles[i].COLOR.g =(i + 10 / 256) % 256;
        Particles[i].COLOR.b =(i + 10 / 65536) % 256;
    }
    return Particles;
}


AlphaParticle *simulate_particles(const Nucleus nuclei[], AlphaParticle *Particles, int particle_count, int nucleus_count, double dt, int max_steps) {
    int i;
    if (Particles== NULL || particle_count <= 0) {
        return NULL;
    }
    for (i=0;i<particle_count;i++) {
    numerical_calculations(&Particles[i], max_steps, nuclei, nucleus_count, dt);
    }
    return Particles;
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
    Nucleus *nuclei = NULL;
    AlphaParticle *Particles; // Array that will store Particles data
    int particle_count = 0;
    double dt;
    int max_steps, WINDOW_HEIGHT, WINDOW_WIDTH;

    // Phase 2: User Input & File Handling
    read_parameters();       // Zsolt
    load_nuclei();           // Ábel

    // Phase 1: Physics & Simulation
    initialize_particles(particle_count, &Particles->position, &Particles->velocity, WINDOW_HEIGHT);  // Pityesz
    simulate_particles(nuclei, &Particles, particle_count, nucleus_count, dt, max_steps);    // Pityesz
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