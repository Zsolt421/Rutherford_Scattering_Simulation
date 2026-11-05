#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*
extensions: Colorful Comments, indent-rainbow,  Party Owl '84, vscode-pets
! Note: Image dimensions are measured in pixels (px), and
! atomic number (Z) is dimensionless.

 ? Recommended simulation parameters (SI units)
 ?
 ? Parameter          Value
 ? -----------------------------
 ? Image width        800 pixels
 ? Image height       800 pixels
 ? Physical width     1e-9 m
 ? Physical height    1e-9 m
 ? Starting X         0 m
 ? Initial velocity   1e7 m/s
 ? Time step (dt)     1e-19 s
 ? Maximum steps      1000
 */


/* ------------------------- 
* Definition of structures     
 ------------------------- */

/* //can be used to store position(x, y), velocity(vx, vy), force(Fx, Fy), accelaration(ax,ay) */
typedef struct
{
    double x;
    double y;
} Vector2;
/* //position, atomic number (positive charge, defaults to 79 if missing) */
typedef struct {
    double x;
    double y;
    int Z;
} Nucleus;
/* //color definition for unique colors */
typedef struct {
    unsigned char r;
    unsigned char g;
    unsigned char b;
} RGB;
/* //stores data for each alpha particle */
typedef struct {
    Vector2 position;
    Vector2 velocity;
    RGB COLOR;
} AlphaParticle;

typedef struct {
    int particle_count;
    int image_width;
    int image_height;
    double physical_width;
    double physical_height;
    double starting_x;
    double initial_velocity;
    double dt;
    int max_steps;
    char nuclei_filename[256];
    char output_filename[256];
} SimulationParameters;

/* ------------------------- 
* Definition of constants      
 ------------------------- */
#define PI 3.14159265358979
#define EPSILON_0 8.8541878188e-12
#define COULOMB_CONSTANT (1.0 / (4.0 * PI * EPSILON_0))

#define ELEMENTARY_CHARGE 1.602176634e-19
#define ALPHA_CHARGE (2.0 * ELEMENTARY_CHARGE)
#define ALPHA_MASS 6.644657e-27


Vector2 vec_add(Vector2 a, Vector2 b) {
    Vector2 res;
    res.x = a.x + b.x;
    res.y = a.y + b.y;
    return res;
}

Vector2 vec_scale(Vector2 v, double factor) {
    Vector2 res;
    res.x = v.x * factor;
    res.y = v.y * factor;
    return res;
}

/* ------------------------- 
* Functions of phase 3       
 ------------------------- */
/*
Calculates the Coulomb force between one alpha particle and one nucleus.
The position is the position of an Alpha particle.
Returns the Coulomb force as a Vector2 (Fx, Fy), measured in newtons.
 */
Vector2 calculate_coulomb_force(Vector2 position, Nucleus nucleus){
    double Fx, Fy, dx, dy, d;
    Vector2 v;

    dx = position.x - nucleus.x;
    dy = position.y - nucleus.y;
    d = sqrt(dx * dx + dy * dy);
    /* the coulomb force multiplied by the direction components, sin(dy/d) or the cos(dx/d), but bc we dont know the angle we use the dx formula */
    Fx = COULOMB_CONSTANT * (ALPHA_CHARGE * nucleus.Z * ELEMENTARY_CHARGE) * dx / (d * d * d);
    Fy = COULOMB_CONSTANT * (ALPHA_CHARGE * nucleus.Z * ELEMENTARY_CHARGE) * dy / (d * d * d);

    v.x = Fx;
    v.y = Fy;

    return v;
}
/* 
todo: function cant handle very small distances btw particles becouse very small nonzero d can lead to division by 0 
*/

Vector2 calculate_total_force(Vector2 position, const Nucleus nuclei[], int nucleus_count) {
    /* F(r) = sum of Fi(r) */
    Vector2 Sum_of_forces = {0,0};
    Vector2 f;
    int i;

    for (i = 0; i < nucleus_count; i++) {
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
    Vector2 k1_r, k1_v, F1;
    Vector2 k2_r, k2_v, pos2, F2;
    Vector2 k3_r, k3_v, pos3, F3;
    Vector2 k4_r, k4_v, pos4, F4;
    Vector2 sum_r, sum_v, r_new, v_new;
    
    while (t <= max_steps){
        
        
        /* first step */
        k1_r = (*alpha).velocity;  /* k1_r = v */
        
        F1 = calculate_total_force((*alpha).position, nuclei, nucleus_count);
        k1_v = vec_scale(F1, 1.0/ ALPHA_MASS); /* k1_v = F(r)/m */


        /* second step */
        k2_r = vec_add((*alpha).velocity, vec_scale(k1_v, 0.5*dt)); /* k2_r = v + 0.5 * dt * k1_v*/

        pos2 = vec_add((*alpha).position, vec_scale(k1_r, 0.5*dt)); /* r + 0.5 * dt * k1_r */
        F2 = calculate_total_force(pos2, nuclei, nucleus_count);       /* F(r + 0.5 * dt * k1_r) */
        k2_v = vec_scale(F2, 1.0/ ALPHA_MASS); /* Y k2_v = F(r + 0.5 * dt * k1_r) / m */


        /* third step */
        k3_r = vec_add((*alpha).velocity, vec_scale(k2_v, 0.5*dt)); /* k3_r = v + 0.5 * dt * k2_v*/

        pos3 = vec_add((*alpha).position, vec_scale(k2_r, 0.5*dt)); /* r + 0.5 * dt * k2_r */
        F3 = calculate_total_force(pos3, nuclei, nucleus_count);       /* F(r + 0.5 * dt * k2_r) */
        k3_v = vec_scale(F3, 1.0/ ALPHA_MASS); /* k3_v = F(r + 0.5 * dt * k2_r) / m */


        /* fourth step */
        k4_r = vec_add((*alpha).velocity, vec_scale(k3_v, dt)); /* k4_r = v + dt * k3_v*/

        pos4 = vec_add((*alpha).position, vec_scale(k3_r, dt)); /* r + dt * k3_r */
        F4 = calculate_total_force(pos4, nuclei, nucleus_count);       /* F(r + dt * k3_r) */
        k4_v = vec_scale(F4, 1.0 / ALPHA_MASS); /* k4_v = F(r + dt * k3_r) / m */


        /* Final update */
        sum_r.x = (k1_r.x + 2*k2_r.x + 2*k3_r.x + k4_r.x);
        sum_r.y = (k1_r.y + 2*k2_r.y + 2*k3_r.y + k4_r.y);

        r_new = vec_add((*alpha).position,vec_scale(sum_r, dt/6)); /*r_new = r + dt * (k1_r + 2*k2_r + 2*k3_r + k4_r) / 6 */
        
        sum_v.x = (k1_v.x + 2*k2_v.x + 2*k3_v.x + k4_v.x);
        sum_v.y = (k1_v.y + 2*k2_v.y + 2*k3_v.y + k4_v.y);

        v_new = vec_add((*alpha).velocity,vec_scale(sum_v, dt/6)); /*v_new = v + dt * (k1_v + 2*k2_v + 2*k3_v + k4_v) / 6 */

        /* SAVE IT INTO A FILE (PPM?) */

        draw_trajectories((*alpha).position, r_new);

        /* increase step */
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
        Particles[i].velocity.x = velocity[i].x; /* Zsolt fell free to change to name of this but pls make sure */
        Particles[i].velocity.y = 0.0;            /* to pay attention that here I assume initial_velocity in a way such that it has only x component */
        Particles[i].position.y = (i + 0.5) * WINDOW_HEIGHT / particle_count; /* We must decide on the height of the playing field (here I used it as double) */
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


/* ------------------------- 
* Functions of phase 2       
 ------------------------- */

SimulationParameters read_parameters(void) {
    SimulationParameters SIMPAR; /* SIMULATION PARAMETERS */

    printf("Enter image width (pixels): ");
    scanf("%d", &SIMPAR.image_width);

    printf("Enter image height (pixels): ");
    scanf("%d", &SIMPAR.image_height);

    printf("Enter physical width (m)(optimal around 10^-9): ");
    scanf("%lf", &SIMPAR.physical_width);

    printf("Enter physical height (m)(optimal around 10^-9): ");
    scanf("%lf", &SIMPAR.physical_height);

    printf("Enter number of alpha particles: ");
    scanf("%d", &SIMPAR.particle_count);

    printf("Enter starting X position (m): ");
    scanf("%lf", &SIMPAR.starting_x);

    printf("Enter initial velocity (m/s): ");
    scanf("%lf", &SIMPAR.initial_velocity);

    printf("Enter time step dt (s): ");
    scanf("%lf", &SIMPAR.dt);

    printf("Enter maximum number of steps: ");
    scanf("%d", &SIMPAR.max_steps);

    printf("Enter nuclei input filename: ");
    scanf("%255s", SIMPAR.nuclei_filename);

    printf("Enter output PPM filename: ");
    scanf("%255s", SIMPAR.output_filename);

    return SIMPAR;
}


int validate_parameters(SimulationParameters SP) {
    if (SP.image_height > 0 && SP.image_width > 0 && SP.physical_height > 0 && SP.physical_width > 0) {
        if (SP.particle_count > 0 && SP.starting_x > 0 && SP.starting_x < 1000 && SP.initial_velocity > 0 && SP.dt > 0 && SP.dt < 100 && SP.max_steps > 0) {
            //!starting_x
            return 1;
        }
    }
    return 0;
}


int validate_file(SimulationParameters SP) {
    FILE *file;

    if (SP.nuclei_filename[0] == '\0' || SP.output_filename[0] == '\0') {
        printf("Empty filename.\n");
        return 1;
    }

    if (strcmp(SP.nuclei_filename, SP.output_filename) == 0) {
        printf("Input and output filenames are identical.\n");
        return 1;
    }

    file = fopen(SP.nuclei_filename, "r");

    if (file == NULL) {
        printf("Cannot open nuclei file.\n");
        return 1;
    }

    fclose(file);

    return 0;
}


Nucleus *read_nuclei(const char *filename, int *out_count) {
    FILE *file;
    Nucleus *p = NULL;
    Nucleus *temp;
    char line[256];
    int count = 0;
    double x, y;
    int Z;
    int read_count;

    file = fopen(filename, "r");
    if (file == NULL) {
        *out_count = 0;
        return NULL;
    }

    while (fgets(line, sizeof(line), file) != NULL) {

        read_count = sscanf(line, "%lf %lf %d", &x, &y, &Z);

        if (read_count >= 2) {
            temp = (Nucleus *)realloc(p, sizeof(Nucleus) * (count + 1));
            if (temp == NULL) {
                free(p);
                fclose(file);
                *out_count = 0;
                return NULL;
            }
            p = temp;

            p[count].x = x;
            p[count].y = y;
            p[count].Z = (read_count == 3) ? Z : 79;

            count++;
        }
    }

    fclose(file);
    *out_count = count;
    return p;
}

/* ------------------------- 
* Functions of phase 3       
 ------------------------- */

int physical_to_pixel_x(double x, double physical_width, int image_width) {
    int px = round(x / physical_width * (image_width - 1));
    return px;
} //& x is the position of the object


int pixel_to_physical(double px, double physical_width, int image_width) {
    int x = px / physical_width * (image_width - 1);
    return x;
} //& px is the pixel position

/* ------------------------- 
* Main                      
 ------------------------- */
int main() {
    int nucleus_count = 0; /* will be scanned from file, 0 is only placeholder */
    Nucleus *nuclei = NULL;
    AlphaParticle *Particles; /* Array that will store Particles data */

    /* Phase 2: User Input & File Handling */
    SimulationParameters SIM_PAR = read_parameters();     /* Zsolt */
    if (1 == validate_parameters(SIM_PAR)){
        printf("Simulation parameters validated.\n");
    }
    else{
        fprintf(stderr,"There is one or more invalid parameters.\n");
        return 1;
    }
    if (validate_file == 1){
        return 1;
    }

    nuclei = read_nuclei(SIM_PAR.nuclei_filename, &nucleus_count);         /* Ábel - nuclei[0] -> nuclei[count-1] */
    if (nuclei == NULL && nucleus_count == 0) {
        fprintf(stderr,"something went wrong with the file.\n");
        return 1;
    }

    /* Phase 1: Physics & Simulation */
    Particles = initialize_particles(SIM_PAR.particle_count, &Particles->position, &Particles->velocity, SIM_PAR.image_height); /* Pityesz */
    if (Particles == NULL) {
        fprintf(stderr,"something went wrong with with initialising the particles\n");
        return 1;
    }
    simulate_particles(nuclei, &Particles, SIM_PAR.particle_count, nucleus_count, SIM_PAR.dt, SIM_PAR.max_steps); /* Pityesz */
    /* Inside simulate_particles(), Pityesz will call Ábel's rk4_step(), which will use Zsolt's three physics functions. */

    /* Phase 3: Visualization & Output */
    create_image();        /* Zsolt */
    draw_nuclei();         /* Ábel */
    draw_trajectories();   /* Ábel */
    save_ppm();            /* Pityesz */

    /* Clean up */
    free(nuclei);
    free_memory();         /* Shared responsibility */

    return 0;
}