#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <SDL/SDL.h>
#include <SDL/SDL_main.h>


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

typedef struct {
    int x;
    int y
} Pixel2;

/* -------------------------
* Definition of constants
 ------------------------- */
#define PI 3.14159265358979
#define EPSILON_0 8.8541878188e-12
#define COULOMB_CONSTANT (1.0 / (4.0 * PI * EPSILON_0))

#define ELEMENTARY_CHARGE 1.602176634e-19
#define ALPHA_CHARGE (2.0 * ELEMENTARY_CHARGE)
#define ALPHA_MASS 6.644657e-27

/* -------------------------
* Functions
 ------------------------- */

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


int physical_to_pixel_x(double x, double physical_width, int image_width) {
    int px = round(x / physical_width * (image_width - 1));
    return px;
} //x is the position of the object


int physical_to_pixel_y(double y, double physical_height, int image_height) {
    int py = round(y / physical_height * (image_height - 1));
    return py;
} //x is the position of the object


/* -------------------------
* Functions of phase 3
 ------------------------- */


Vector2 calculate_coulomb_force(Vector2 position, Nucleus nucleus){
    long double Fx, Fy, dx, dy, d;
    Vector2 v;

    dx = position.x - nucleus.x;
    dy = position.y - nucleus.y;
    d = sqrt(dx * dx + dy * dy);
    /* the coulomb force multiplied by the direction components, sin(dy/d) or the cos(dx/d), but bc we dont know the angle we use the dx formula */
    Fx = COULOMB_CONSTANT * (ALPHA_CHARGE * nucleus.Z * ELEMENTARY_CHARGE) * dx / (d * d * d);
    Fy = COULOMB_CONSTANT * (ALPHA_CHARGE * nucleus.Z * ELEMENTARY_CHARGE) * dy / (d * d * d);

    v.x = (double)Fx;
    v.y = (double)Fy;

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

        r_new = vec_add((*alpha).position, vec_scale(sum_r, dt/6)); /*r_new = r + dt * (k1_r + 2*k2_r + 2*k3_r + k4_r) / 6 */

        sum_v.x = (k1_v.x + 2*k2_v.x + 2*k3_v.x + k4_v.x);
        sum_v.y = (k1_v.y + 2*k2_v.y + 2*k3_v.y + k4_v.y);

        v_new = vec_add((*alpha).velocity, vec_scale(sum_v, dt/6)); /*v_new = v + dt * (k1_v + 2*k2_v + 2*k3_v + k4_v) / 6 */

        /* SAVE IT INTO A FILE (PPM?) */

        draw_trajectories((*alpha).position, r_new);
        //draw_trajectories(Vector2 prev_pos, Vector2 pos, RGB color, int physical_width, int image_width, int physical_height, int image_height, RGB *image);

        /* increase step */
        (*alpha).position = r_new;
        (*alpha).velocity = v_new;
        t++;
    }

    return;
}


AlphaParticle *initialize_particles(SimulationParameters SIMPAR) {
    AlphaParticle *Particles;
    int i;
    if (SIMPAR.particle_count <= 0 || SIMPAR.initial_velocity == NULL) {
        return NULL;
        }
    Particles = (AlphaParticle *)malloc(sizeof(AlphaParticle) * SIMPAR.particle_count);
    if (Particles == NULL) {
        return NULL;
    }
    for (i = 0; i < SIMPAR.particle_count; i++) {
        Particles[i].position.x = SIMPAR.starting_x;
        Particles[i].velocity.x = SIMPAR.initial_velocity[i]; /* Zsolt fell free to change to name of this but pls make sure */
        Particles[i].velocity.y = 0.0;            /* to pay attention that here I assume initial_velocity in a way such that it has only x component */
        Particles[i].position.y = ((i + 0.5) * SIMPAR.physical_height / SIMPAR.particle_count); /* We must decide on the height of the playing field (here I used it as double) */
        Particles[i].COLOR.r = (i + 10) % 256;
        Particles[i].COLOR.g = ((i + 10) / 256) % 256;
        Particles[i].COLOR.b = ((i + 10) / 65536) % 256;
    }
    return Particles;
}


AlphaParticle *simulate_particles(const Nucleus nuclei[], AlphaParticle *Particles, int nucleus_count, SimulationParameters SIMPAR) {
    int i;
    if (Particles== NULL || SIMPAR.particle_count <= 0 || SIMPAR.dt  <= 0.0 || SIMPAR.max_steps <= 0 || nucleus_count < 0 || (nucleus_count > 0 && nuclei == NULL)) {
        return NULL;
    }
    for (i = 0; i < SIMPAR.particle_count; i++) {
    numerical_calculations(&Particles[i], SIMPAR.max_steps, nuclei, nucleus_count, SIMPAR.dt);
    }
    return Particles;
}


/* -------------------------
* Functions of phase 2
 ------------------------- */


SimulationParameters read_parameters(void) {              //Here I can't use return 1 to exit because my return tyoe uns SMP so I use exit(EXIT_FAILURE) wich kills of the program
    SimulationParameters SIMPAR = {0}; /* SIMULATION PARAMETERS */

    printf("Enter image width (pixels): ");
    if (scanf("%d", &SIMPAR.image_width) !=1) {
        fprintf(stderr, "Invalid image width\n");
        exit(EXIT_FAILURE);
    }

    printf("Enter image height (pixels): ");
    if (scanf("%d", &SIMPAR.image_height) != 1) {
        fprintf(stderr, "Invalid image height\n");
        exit(EXIT_FAILURE);
    }

    printf("Enter physical width (m)(optimal around 10^-9): ");
    if (scanf("%lf", &SIMPAR.physical_width) != 1) {
        fprintf(stderr, "Invalid physical width\n");
        exit(EXIT_FAILURE);
    }

    printf("Enter physical height (m)(optimal around 10^-9): ");
    if (scanf("%lf", &SIMPAR.physical_height) != 1) {
        fprintf(stderr, "Invalid physical height\n");
        exit(EXIT_FAILURE);
    }

    printf("Enter number of alpha particles: "); //Error if invalid particle count
    if (scanf("%d", &SIMPAR.particle_count)!=1 || SIMPAR.particle_count <= 0) {
        fprintf(stderr, "Invalid number of alpha particles\n");
        exit(EXIT_FAILURE);
    }

    printf("Enter starting X position (m): ");
    if (scanf("%lf", &SIMPAR.starting_x) != 1) {
        fprintf(stderr, "Invalid starting X position\n");
        exit(EXIT_FAILURE);
    }

    int i;  //Creating velocity array
    SIMPAR.initial_velocity = malloc(SIMPAR.particle_count * sizeof(double));
    if (SIMPAR.initial_velocity == NULL) {
        fprintf(stderr, "Error allocating initial velocity\n"); //Checking if user entered 0 initial
        exit(EXIT_FAILURE);
    }
    for (i=0;i<SIMPAR.particle_count;i++) {
    printf("Enter initial velocity %d (m/s): ",i+1); //allocating velocities
    if (scanf("%lf", &SIMPAR.initial_velocity[i]) != 1) { //Checking if entered velocity is a correct type
        fprintf(stderr, "Invalid initial velocity\n");
        free(SIMPAR.initial_velocity);
        exit(EXIT_FAILURE);
        }
    }
    printf("Enter time step dt (s): ");
    if (scanf("%lf", &SIMPAR.dt) != 1) {
        fprintf(stderr, "Invalid time step input\n");
        free(SIMPAR.initial_velocity);
        exit(EXIT_FAILURE);
    }

    printf("Enter maximum number of steps: ");
    if (scanf("%d", &SIMPAR.max_steps) != 1) {
        fprintf(stderr, "Invalid number of steps\n");
        free(SIMPAR.initial_velocity);
        exit(EXIT_FAILURE);
    }

    printf("Enter nuclei input filename: ");
    if (scanf("%255s", SIMPAR.nuclei_filename) != 1) {
        fprintf(stderr, "Invalid filename\n");
        free(SIMPAR.initial_velocity);
        exit(EXIT_FAILURE);
    }

    printf("Enter output PPM filename: ");
    if (scanf("%255s", SIMPAR.output_filename) != 1) {
        fprintf(stderr, "Invalid filename\n");
        free(SIMPAR.initial_velocity);
        exit(EXIT_FAILURE);
    }

    return SIMPAR;
}


int validate_parameters(SimulationParameters SIMPAR) { //Here I use isfinite tha checks if my number is finite only returns 0 if isfinite(INFINITY) ors isfinte(NAN) !isfinite is its inverse so it returns 0 when isfinite would return 1
    int i;

    if (SIMPAR.image_width <= 0 || SIMPAR.image_height <= 0) {
        fprintf(stderr, "Invalid image dimensions.\n");
        return 1;
    }

    if (!isfinite(SIMPAR.physical_width) || !isfinite(SIMPAR.physical_height) || SIMPAR.physical_width <= 0.0 ||
        SIMPAR.physical_height <= 0.0) {
        fprintf(stderr, "Invalid physical dimensions.\n");
        return 1;
        }

    if (SIMPAR.particle_count <= 0 || SIMPAR.initial_velocity == NULL) {
        fprintf(stderr, "Invalid particle count or velocity array.\n");
        return 1;
        }

    if (!isfinite(SIMPAR.starting_x) || SIMPAR.starting_x < 0.0 || SIMPAR.starting_x >= SIMPAR.physical_width) {
        fprintf(stderr, "Invalid starting X position.\n");
        return 1;
        }

    for (i = 0; i < SIMPAR.particle_count; i++) {
        if (!isfinite(SIMPAR.initial_velocity[i]) || SIMPAR.initial_velocity[i] <= 0.0) {
            fprintf(stderr,"Invalid initial velocity for particle %d.\n",i + 1);
            return 1;
            }
    }

    if (!isfinite(SIMPAR.dt) || SIMPAR.dt <= 0.0 || SIMPAR.dt >= 1.0) {
        fprintf(stderr, "Invalid time step.\n");
        return 1;
        }

    if (SIMPAR.max_steps <= 0) {
        fprintf(stderr, "Invalid maximum number of steps.\n");
        return 1;
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

        if (read_count < 2 || !isfinite(x) || !isfinite(y)) {
            fprintf(stderr,"Invalid nuclei data on line %d.\n",count+1);
            free(p);
            fclose(file);
            *out_count = 0;
            return NULL;
        }
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

        if (p[count].Z<=0) {
            fprintf(stderr,"Invalid nuclei data on line %d.\n",count+1);
            free(p);
            fclose(file);
            *out_count = 0;
            return NULL;
        }
        count++;
    }

    fclose(file);
    *out_count = count;
    return p;
}


/* -------------------------
* Functions of phase 3
 ------------------------- */


RGB *create_image(int width, int height) {
    RGB *image;
    int pixel_count, i;
    pixel_count = (size_t)width * (size_t)height;
    image = malloc(pixel_count * sizeof(RGB));


    if (image == NULL) {
        return NULL;
    }
    for (i = 0; i < pixel_count; i++) {
        image[i].r = 0;
        image[i].g = 0;
        image[i].b = 0;
    }

    return image;
}


int color_pixel(int image_width, int iy, int ix, RGB *image, RGB color) {
    int position = round(iy * image_width + ix);
    image[position].r = color.r;
    image[position].g = color.g;
    image[position].b = color.b;
    return 1;
}

//!not functional
int draw_trajectories(Vector2 prev_pos, Vector2 pos, RGB color, int physical_width, int image_width, int physical_height, int image_height, RGB *image) {
    int ix, iy;
    Pixel2 pp;
    Pixel2 p;
    pp.x = physical_to_pixel_x(prev_pos.x, physical_width, image_width);
    p.x = physical_to_pixel_x(pos.x, physical_width, image_width);
    pp.y = physical_to_pixel_y(prev_pos.y, physical_height, image_height);
    p.y = physical_to_pixel_y(pos.y, physical_height, image_height);

    int dx = p.x - pp.x;
    int dy = p.y - pp.y;

    int steps;

    if (abs(dx) > abs(dy)) {
        steps = abs(dx);
    }
    else {
        steps = abs(dy);
    }

    double x = pp.x;
    double y = pp.y;

    if (steps == 0) {
        color_pixel(image_width, pp.y, pp.x, image, color);
        return 1;
    }

    double x_step = (double)dx / steps;
    double y_step = (double)dy / steps;

    for (int i = 0; i <= steps; i++) {
        ix = (int)round(x);
        iy = (int)round(y);

        if (ix >= 0 && ix < image_width && iy >= 0 && iy < image_height) {
            color_pixel(image_width, iy, ix, image, color);
        }

        x += x_step;
        y += y_step;
    }

    return 1;

}

int save_ppm(const char *filename, const RGB *image, int width, int height) {
    FILE *ppm=fopen(filename, "w");
    if (ppm == NULL) {
        fprintf(stderr,"Unable to open file for writing.\n");
        return 1;
    }
    fprintf(ppm,"P3\n %d %d\n 255\n", width, height);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            fprintf(ppm, "%d %d %d ", image[i*width+j].r,image[i*width+j].g,image[i*width+j].b);
        }
        fprintf(ppm, "\n");
    }
    fclose(ppm);
    return 0;
}


/* -------------------------
* Main
 ------------------------- */


int main(int argc, char *argv[]) {
    int nucleus_count = 0; /* will be scanned from file, 0 is only placeholder */
    Nucleus *nuclei = NULL;
    AlphaParticle *Particles = NULL; /* Array that will store Particles data */

    /* Phase 2: User Input & File Handling */
    SimulationParameters SIM_PAR = read_parameters();     /* Zsolt */
    if (1 == validate_parameters(SIM_PAR)){
        printf("Simulation parameters validated.\n");
    }
    else{
        fprintf(stderr,"There is one or more invalid parameters.\n");
        free(SIM_PAR.initial_velocity);
        return 1;
    }

    if (validate_file(SIM_PAR) == 1) {
        free(SIM_PAR.initial_velocity);
        return 1;
    }

    RGB *image = create_image(SIM_PAR.image_width, SIM_PAR.image_height);  /* Zsolt */

    nuclei = read_nuclei(SIM_PAR.nuclei_filename, &nucleus_count);         /* Ábel - nuclei[0] -> nuclei[count-1] */
    if (nuclei == NULL && nucleus_count == 0) {
        fprintf(stderr,"something went wrong with the file.\n");
        return 1;
    }

    /* Phase 1: Physics & Simulation */
    Particles = initialize_particles(SIM_PAR); /* Pityesz */
    if (Particles == NULL) {
        fprintf(stderr,"something went wrong with with initialising the particles\n");
        return 1;
    }
    simulate_particles(nuclei, Particles, nucleus_count, SIM_PAR); /* Pityesz */
    /* Inside simulate_particles(), Pityesz will call Ábel's rk4_step(), which will use Zsolt's three physics functions. */

    /* Phase 3: Visualization & Output */
    draw_nuclei();         /* Ábel */
    save_ppm();            /* Pityesz */

    /* Clean up */
    free(nuclei);/* Shared responsibility */

    return 0;
}