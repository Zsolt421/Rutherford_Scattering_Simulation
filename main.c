#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <SDL/SDL.h>
#undef main
#include <SDL/SDL_main.h>

#include "simulation.h"


/* -------------------------
* Functions
 ------------------------- */


/* Component-wise vector addition, shared by the integration calculations. */
Vector2 vec_add(Vector2 a, Vector2 b) {
    Vector2 res;
    res.x = a.x + b.x;
    res.y = a.y + b.y;
    return res;
}


/* Multiply both vector components by the same scalar. */
Vector2 vec_scale(Vector2 v, double factor) {
    Vector2 res;
    res.x = v.x * factor;
    res.y = v.y * factor;
    return res;
}


/* Map a physical coordinate in the simulation region to the nearest pixel.
The region's endpoints map to pixel 0 and image_width - 1 respectively.
 */
int physical_to_pixel_x(double x, double physical_width, int image_width) {
    int px = round(x / physical_width * (image_width - 1));
    return px;
}


int physical_to_pixel_y(double y, double physical_height, int image_height) {
    int py = round(y / physical_height * (image_height - 1));
    return py;
}

//Progress bar
void show_progress(long long completed, long long total) {
    int width = 40;
    int filled = (int)((double)completed / total * width);
    int percent = (int)((double)completed / total * 100);

    printf("\rSimulating [");

    for (int i = 0; i < width; i++) {
        if (i < filled) printf("=");
        else if (i == filled && completed < total) printf(">");
        else printf(" ");
    }

    printf("] %3d%%", percent);
    fflush(stdout);
}


/* -------------------------
 * Functions of phase 2: Input and validation
 ------------------------- */

/* Append .txt to a filename if the extension is missing. */
char *txt_adder(char *filename) {
    size_t len = strlen(filename);

    /* Keep the original name if it already ends with .txt. */
    if (len < 4 || strcmp(filename + len - 4, ".txt") != 0) {
        if (len > 251) {
            fprintf(stderr, "Filename is too long.\n");
            return NULL;
        }
        strcat(filename, ".txt");
    }
    return filename;
}

/* Append .ppm to a filename if the extension is missing. */
char *ppm_adder(char *filename) {
    size_t len = strlen(filename);

    /* Keep the original name if it already ends with .ppm. */
    if (len < 4 || strcmp(filename + len - 4, ".ppm") != 0) {
        if (len > 251) {
            fprintf(stderr, "Filename is too long.\n");
            return NULL;
        }
        strcat(filename, ".ppm");
    }
    return filename;
}


/* Read the interactive configuration. File extensions are appended when
omitted; parameter ranges and file checks are handled by the validation
and input-reading routines.
 */
SimulationParameters read_parameters_demo(void) {
    SimulationParameters SIMPAR = {0};
    SIMPAR.image_width = 1200;
    SIMPAR.image_height = 1200;
    SIMPAR.physical_width = 1e-12;
    SIMPAR.physical_height = 1e-12;
    SIMPAR.particle_count = 50;
    SIMPAR.starting_x = 0;
    SIMPAR.initial_velocity_min = 5e6;
    SIMPAR.initial_velocity_max = 8e7;
    SIMPAR.dt = 1e-22;
    SIMPAR.max_steps = 100000;
    strcpy(SIMPAR.output_filename, "output.ppm");
    printf("Demo parameters loaded\n");
    return SIMPAR;
}


Validate read_parameters(SimulationParameters *SIMPAR) {
    *SIMPAR = (SimulationParameters){0}; /* Start with empty names and zeroed fields. */

    printf("Enter image width (pixels): ");
    fflush(stdout);
    if (scanf("%d", &SIMPAR->image_width) !=1) {
        fprintf(stderr, "Invalid image width\n");
        return INVALID;
    }

    printf("Enter image height (pixels): ");
    fflush(stdout);
    if (scanf("%d", &SIMPAR->image_height) != 1) {
        fprintf(stderr, "Invalid image height\n");
        return INVALID;
    }

    printf("Enter physical width (m)(optimal around 10^-9): ");
    fflush(stdout);
    if (scanf("%lf", &SIMPAR->physical_width) != 1) {
        fprintf(stderr, "Invalid physical width\n");
        return INVALID;
    }

    printf("Enter physical height (m)(optimal around 10^-9): ");
    fflush(stdout);
    if (scanf("%lf", &SIMPAR->physical_height) != 1) {
        fprintf(stderr, "Invalid physical height\n");
        return INVALID;
    }

    printf("Enter number of alpha particles: ");
    fflush(stdout);
    if (scanf("%d", &SIMPAR->particle_count)!=1 || SIMPAR->particle_count <= 0) {
        fprintf(stderr, "Invalid number of alpha particles\n");
        return INVALID;
    }

    printf("Enter starting X position (m): ");
    fflush(stdout);
    if (scanf("%lf", &SIMPAR->starting_x) != 1) {
        fprintf(stderr, "Invalid starting X position\n");
        return INVALID;
    }

    printf("Enter velocity range minimum (m/s): ");
    fflush(stdout);
    if (scanf("%lf", &SIMPAR->initial_velocity_min) != 1) {
        fprintf(stderr, "Invalid velocity range minimum\n");
        return INVALID;
    }

    printf("Enter velocity range maximum (m/s): ");
    fflush(stdout);
    if (scanf("%lf", &SIMPAR->initial_velocity_max) != 1) {
        fprintf(stderr, "Invalid velocity range maximum\n");
        return INVALID;
    }

    printf("Enter time step dt (s): ");
    fflush(stdout);
    if (scanf("%lf", &SIMPAR->dt) != 1) {
        fprintf(stderr, "Invalid time step input\n");
        return INVALID;
    }

    printf("Enter maximum number of steps: ");
    fflush(stdout);
    if (scanf("%d", &SIMPAR->max_steps) != 1) {
        fprintf(stderr, "Invalid number of steps\n");
        return  INVALID;
    }

    printf("Enter nuclei input filename (with or without .txt): ");
    fflush(stdout);
    if (scanf("%251s", SIMPAR->nuclei_filename) != 1) {
        fprintf(stderr, "Invalid filename scan\n");
        return INVALID;
    }

    printf("Enter output PPM filename (with or without .ppm): ");
    fflush(stdout);
    if (scanf("%251s", SIMPAR->output_filename) != 1) {
        fprintf(stderr, "Invalid filename\n");
        return INVALID;
    }

    /* Use the same filename helpers as in the file input mode. */
    if (txt_adder(SIMPAR->nuclei_filename) == NULL) {
        return INVALID;
    }
    if (ppm_adder(SIMPAR->output_filename) == NULL) {
        return INVALID;
    }

    return VALID;
}


Validate read_parameters_from_file(SimulationParameters *SIMPAR) {
    *SIMPAR = (SimulationParameters){0}; /* Start with empty names and zeroed fields. */
    char input_file_name_alpha[256], output_file_name[256], input_file_name_nuclei[256];

    printf("Enter input file name of simulation parameters (with or without .txt) (Max 255 characters): ");
    fflush(stdout);
    if (scanf("%255s", input_file_name_alpha) != 1) {
        fprintf(stderr, "Invalid input file name of simulation parameters.\n");
        return INVALID;
    }
    if (input_file_name_alpha[0] == '\0') {
        fprintf(stderr, "Input file name of simulation parameters cannot be empty.\n");
        return INVALID;
    }
    /* Append .txt if missing, without exceeding the filename array size. */
    if (txt_adder(input_file_name_alpha) == NULL) {
        return INVALID;
    }

    printf("Enter input file name of nuclei (with or without .txt) (Max 255 characters): ");
    fflush(stdout);
    if (scanf("%255s", input_file_name_nuclei) != 1) {
        fprintf(stderr, "Invalid input file name of nuclei.\n");
        return INVALID;
    }
    if (input_file_name_nuclei[0] == '\0') {
        fprintf(stderr, "Input file name of nuclei cannot be empty.\n");
        return INVALID;
    }
    /* Append .txt if missing, without exceeding the filename array size. */
    if (txt_adder(input_file_name_nuclei) == NULL) {
        return INVALID;
    }

    /* Store the nuclei filename for read_nuclei(). */
    /* Include the terminating null character when copying the name. */
    for (size_t i = 0; i <= strlen(input_file_name_nuclei); i++) {
        SIMPAR->nuclei_filename[i] = input_file_name_nuclei[i];
    }

    printf("Enter output file name (with or without .ppm) (Max 255 characters): ");
    fflush(stdout);
    if (scanf("%255s", output_file_name) != 1) {
        fprintf(stderr, "Invalid output file name.\n");
        return INVALID;
    }
    if (output_file_name[0] == '\0') {
        fprintf(stderr, "Output file name cannot be empty.\n");
        return INVALID;
    }
    /* Append .ppm if missing, without exceeding the filename array size. */
    if (ppm_adder(output_file_name) == NULL) {
        return INVALID;
    }

    /* Store the output filename for save_ppm(). */
    /* Include the terminating null character when copying the name. */
    for (size_t i = 0; i <= strlen(output_file_name); i++) {
        SIMPAR->output_filename[i] = output_file_name[i];
    }

    /* Numerical parameters and nuclei positions are stored in separate files. */
    FILE *file = fopen(input_file_name_alpha, "r");
    if (file == NULL) {
        fprintf(stderr, "Could not open %s for reading.\n", input_file_name_alpha);
        return INVALID;
    }

    /* Read the ten numerical simulation parameters from the input file. */
    if (fscanf(file, "%d %d %lf %lf %d %lf %lf %lf %lf %d",
        &SIMPAR->image_width,
        &SIMPAR->image_height,
        &SIMPAR->physical_width,
        &SIMPAR->physical_height,
        &SIMPAR->particle_count,
        &SIMPAR->starting_x,
        &SIMPAR->initial_velocity_min,
        &SIMPAR->initial_velocity_max,
        &SIMPAR->dt,
        &SIMPAR->max_steps) != 10) {
        fprintf(stderr, "Invalid %s format.\n", input_file_name_alpha);
        fclose(file);
        return INVALID;
    }

    fclose(file);
    return VALID;
}


/* Check numeric settings for valid dimensions, ranges, and finite values. */
Validate validate_parameters(SimulationParameters SIMPAR) {

    if (SIMPAR.image_width <= 0 || SIMPAR.image_height <= 0) {
        fprintf(stderr, "Invalid image dimensions.\n");
        return INVALID;
    }

    if (!isfinite(SIMPAR.physical_width) || !isfinite(SIMPAR.physical_height) || SIMPAR.physical_width <= 0.0 ||
        SIMPAR.physical_height <= 0.0) {
        fprintf(stderr, "Invalid physical dimensions.\n");
        return INVALID;
        }

    if (SIMPAR.particle_count <= 0 || SIMPAR.initial_velocity_max < SIMPAR.initial_velocity_min || SIMPAR.initial_velocity_min <= 0 || !isfinite(SIMPAR.initial_velocity_min) || !isfinite(SIMPAR.initial_velocity_max)) {
        fprintf(stderr, "Invalid particle count or velocity range\n");
        return INVALID;
        }

    if (!isfinite(SIMPAR.starting_x) || SIMPAR.starting_x < 0.0 || SIMPAR.starting_x >= SIMPAR.physical_width) {
        fprintf(stderr, "Invalid starting X position.\n");
        return INVALID;
        }

    if (!isfinite(SIMPAR.dt) || SIMPAR.dt <= 0.0 || SIMPAR.dt >= 1.0) {
        fprintf(stderr, "Invalid time step.\n");
        return INVALID;
        }

    if (SIMPAR.max_steps <= 0) {
        fprintf(stderr, "Invalid maximum number of steps.\n");
        return INVALID;
    }

    return VALID;
}


Validate validate_file(SimulationParameters SP) {
    FILE *file;

    /* Reject empty or identical names, then confirm the input file is readable. */
    if (SP.nuclei_filename[0] == '\0' || SP.output_filename[0] == '\0') {
        printf("Empty filename.\n");
        return INVALID;
    }

    if (strcmp(SP.nuclei_filename, SP.output_filename) == 0) {
        printf("Input and output filenames are identical.\n");
        return INVALID;
    }

    file = fopen(SP.nuclei_filename, "r");

    if (file == NULL) {
        printf("Cannot open nuclei file.\n");
        return INVALID;
    }

    fclose(file);

    return VALID;
}


Validate validate_nuclei(SimulationParameters SIMPAR, int read_count, int count, double x, double y, int Z){
    /* Nucleus coordinates must be within the physical region and in its right
    half (x >= physical_width/2). A two-field record uses the default Z.
     */
    if (x < 0.0 || y < 0.0 || x > SIMPAR.physical_width || y > SIMPAR.physical_height || !isfinite(x) || !isfinite(y)) {
        fprintf(stderr, "Nuclei position invalid in %s line %d", SIMPAR.nuclei_filename, count+1);
        return INVALID;
    }

    if (x < SIMPAR.physical_width/2){
        fprintf(stderr, "Nuclei position invalid in %s line %d, nuclei must be on the right side of the screen\n", SIMPAR.nuclei_filename, count+1);
        return INVALID;
    }

    if (read_count == 2){
        return VALID;
    }

    if (Z <= 0 || Z > 92){ /* Atomic number check, exists but not larger than Uranium */
        fprintf(stderr, "Atomic number of the Nuclei invalid in %s line %d", SIMPAR.nuclei_filename, count+1);
        return INVALID;
    }

    return VALID;

}


Nucleus *read_nuclei(SimulationParameters SIMPAR, int *out_count) {
    FILE *file;
    Nucleus *p = NULL;
    Nucleus *temp;
    char line[256];
    int count = 0;
    double x, y;
    int Z;
    int read_count;

    /* Read one whitespace-separated `x y [Z]` record per line. Coordinates
    are in metres; omitted atomic numbers default to 79. Malformed lines
    fail the whole read, while successfully parsed records grow the array.
     */
    file = fopen(SIMPAR.nuclei_filename, "r");
    if (file == NULL) {
        *out_count = 0;
        return NULL;
    }

    while (fgets(line, sizeof(line), file) != NULL) {

        Z = 79;                  /* Default atomic number when the line omits Z. */
        read_count = sscanf(line, "%lf %lf %d", &x, &y, &Z);

        if (read_count >= 2) {

            if (validate_nuclei(SIMPAR,read_count, count,x,y,Z) == INVALID){
                fclose(file);
                free(p);
                *out_count = 0;
                return NULL;
            }

            /* Increase the dynamically allocated array for the next nucleus. */
            temp = (Nucleus *)realloc(p, sizeof(Nucleus) * (count + 1));
            if (temp == NULL) {
                free(p);
                fclose(file);
                *out_count = 0;
                return NULL;
            }
            p = temp;

            p[count].pos.x = x;
            p[count].pos.y = y;
            p[count].Z = (read_count == 3) ? Z : 79;

            count++;
        }

        else {
            fprintf(stderr, "Not enough parameters in %s, line %d", SIMPAR.nuclei_filename, count+1);
            free(p);
            fclose(file);
            *out_count = 0;
            return NULL;
        }
    }

    fclose(file);
    *out_count = count;
    return p;
}


/* -------------------------
 * Functions of phase 1: Physics and simulation
 ------------------------- */

/* Return the electrostatic force on an alpha particle at position from one
positively charged nucleus. The inverse-square Coulomb force is resolved
into x/y components; coincident positions are singular and are not handled
specially here.
 */
Vector2 calculate_coulomb_force(Vector2 position, Nucleus nucleus){
    long double Fx, Fy, dx, dy, d;
    Vector2 v;

    dx = position.x - nucleus.pos.x;
    dy = position.y - nucleus.pos.y;
    d = sqrt(dx * dx + dy * dy);
    /* Resolve the radial force along the displacement vector (dx, dy). */
    Fx = COULOMB_CONSTANT * (ALPHA_CHARGE * nucleus.Z * ELEMENTARY_CHARGE) * dx / (d * d * d);
    Fy = COULOMB_CONSTANT * (ALPHA_CHARGE * nucleus.Z * ELEMENTARY_CHARGE) * dy / (d * d * d);

    v.x = (double)Fx;
    v.y = (double)Fy;

    return v;
}

Vector2 calculate_total_force(Vector2 position, const Nucleus nuclei[], int nucleus_count) {
    /* Superpose the force contributions from all nuclei at the same position. */
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


/* Earlier direct acceleration helper; the integrator below computes F/m at
each Runge-Kutta stage instead.
Vector2 calculate_acceleration(Vector2 position, const Nucleus nuclei[], int nucleus_count) {
    Vector2 acceleration_of_particle = {0,0};
    Vector2 F = calculate_total_force(position, nuclei, nucleus_count);
    acceleration_of_particle.x = F.x / ALPHA_MASS;
    acceleration_of_particle.y = F.y / ALPHA_MASS;
    return acceleration_of_particle;
} */


/* Advance one particle for max_steps using classical fourth-order Runge-Kutta
on dr/dt = v and dv/dt = F(r)/m. Each completed step is drawn before the
particle's position and velocity are updated.
 */
Validate numerical_calculations(AlphaParticle *alpha, const Nucleus nuclei[], int nucleus_count, SimulationParameters SIMPAR, RGB *image, int particle_index){
    int t = 1;
    Vector2 k1_r, k1_v, F1;
    Vector2 k2_r, k2_v, pos2, F2;
    Vector2 k3_r, k3_v, pos3, F3;
    Vector2 k4_r, k4_v, pos4, F4;
    Vector2 sum_r, sum_v, r_new, v_new;

    while (t <= SIMPAR.max_steps) {

        /* First RK4 stage: evaluate the derivatives at the current state. */
        k1_r = (*alpha).velocity;  /** k1_r = v */

        F1 = calculate_total_force((*alpha).position, nuclei, nucleus_count);
        k1_v = vec_scale(F1, 1.0/ ALPHA_MASS); /* k1_v = F(r)/m */


        /* Second stage: estimate the state halfway through the time step. */
        k2_r = vec_add((*alpha).velocity, vec_scale(k1_v, 0.5 * SIMPAR.dt)); /* k2_r = v + 0.5 * dt * k1_v*/

        pos2 = vec_add((*alpha).position, vec_scale(k1_r, 0.5 * SIMPAR.dt)); /* r + 0.5 * dt * k1_r */
        F2 = calculate_total_force(pos2, nuclei, nucleus_count);       /* F(r + 0.5 * dt * k1_r) */
        k2_v = vec_scale(F2, 1.0/ ALPHA_MASS); /* Y k2_v = F(r + 0.5 * dt * k1_r) / m */


        /* Third stage: refine the midpoint estimate using stage two. */
        k3_r = vec_add((*alpha).velocity, vec_scale(k2_v, 0.5 * SIMPAR.dt)); /* k3_r = v + 0.5 * dt * k2_v*/

        pos3 = vec_add((*alpha).position, vec_scale(k2_r, 0.5 * SIMPAR.dt)); /* r + 0.5 * dt * k2_r */
        F3 = calculate_total_force(pos3, nuclei, nucleus_count);       /* F(r + 0.5 * dt * k2_r) */
        k3_v = vec_scale(F3, 1.0/ ALPHA_MASS); /* k3_v = F(r + 0.5 * dt * k2_r) / m */


        /* Fourth stage: estimate the derivatives at the end of the step. */
        k4_r = vec_add((*alpha).velocity, vec_scale(k3_v, SIMPAR.dt)); /* k4_r = v + dt * k3_v*/

        pos4 = vec_add((*alpha).position, vec_scale(k3_r, SIMPAR.dt)); /* r + dt * k3_r */
        F4 = calculate_total_force(pos4, nuclei, nucleus_count);       /* F(r + dt * k3_r) */
        k4_v = vec_scale(F4, 1.0 / ALPHA_MASS); /* k4_v = F(r + dt * k3_r) / m */


        /* Combine the four stages with weights 1:2:2:1 to update the state. */
        sum_r.x = (k1_r.x + 2*k2_r.x + 2*k3_r.x + k4_r.x);
        sum_r.y = (k1_r.y + 2*k2_r.y + 2*k3_r.y + k4_r.y);

        r_new = vec_add((*alpha).position, vec_scale(sum_r, SIMPAR.dt/6)); /**r_new = r + dt * (k1_r + 2*k2_r + 2*k3_r + k4_r) / 6 */

        sum_v.x = (k1_v.x + 2*k2_v.x + 2*k3_v.x + k4_v.x);
        sum_v.y = (k1_v.y + 2*k2_v.y + 2*k3_v.y + k4_v.y);

        v_new = vec_add((*alpha).velocity, vec_scale(sum_v, SIMPAR.dt/6)); /**v_new = v + dt * (k1_v + 2*k2_v + 2*k3_v + k4_v) / 6 */

        if (draw_trajectories((*alpha).position, r_new, (*alpha).COLOR, SIMPAR.physical_width, SIMPAR.image_width, SIMPAR.physical_height, SIMPAR.image_height, image) == INVALID) {
            return INVALID;
        }

        /* Commit the completed step and continue with the new state. */
        (*alpha).position = r_new;
        (*alpha).velocity = v_new;
        t++;
        long long total = (long long)SIMPAR.particle_count * SIMPAR.max_steps;
        long long completed = (long long)particle_index * SIMPAR.max_steps + t - 1;
        long long interval = total / 200;

        if (interval < 1) interval = 1;

        if (completed % interval == 0 || completed == total) {
            show_progress(completed, total);
        }
    }

    return VALID;
}


AlphaParticle *initialize_particles(SimulationParameters SIMPAR) {
    AlphaParticle *Particles;
    int i;
    Particles = malloc(sizeof(AlphaParticle) * SIMPAR.particle_count);
    if (Particles == NULL) {
        return NULL;
    }

    for (i = 0; i < SIMPAR.particle_count; i++) {
        Particles[i].position.x = SIMPAR.starting_x;
        Particles[i].velocity.x = SIMPAR.initial_velocity_min + ((double)rand() / RAND_MAX) * (SIMPAR.initial_velocity_max - SIMPAR.initial_velocity_min);
        /* Particles enter horizontally and are evenly spaced across the field. */
        Particles[i].velocity.y = 0.0;
        Particles[i].position.y = ((i + 0.5) * SIMPAR.physical_height / SIMPAR.particle_count);
        /* Vary each RGB channel by particle index to distinguish trajectories. */
        Particles[i].COLOR.r = 100 + (i * 37LL) % 156;
        Particles[i].COLOR.g = 80 + (i * 67LL) % 176;
        Particles[i].COLOR.b = 100 + (i * 97LL) % 156;
    }
    return Particles;
}


AlphaParticle *simulate_particles(const Nucleus nuclei[], AlphaParticle *Particles, int nucleus_count, SimulationParameters SIMPAR, RGB *image) {
    int i;
    /* Reject inconsistent state before accessing the particle/nucleus arrays. */
    if (Particles== NULL || SIMPAR.particle_count <= 0 || SIMPAR.dt  <= 0.0 || SIMPAR.max_steps <= 0 || nucleus_count < 0 || (nucleus_count > 0 && nuclei == NULL)) {
        return NULL;
    }
    long long total = (long long)SIMPAR.particle_count * SIMPAR.max_steps;
    printf("\x1b[?25l\n");  /* Hide cursor */
    show_progress(0, total);
    for (i = 0; i < SIMPAR.particle_count; i++) {
        if (numerical_calculations(&Particles[i], nuclei, nucleus_count, SIMPAR, image, i) == INVALID) {
            printf("\x1b[?25h\n");
            return NULL;
        }
    }
    show_progress(total, total);
    printf("\x1b[?25h\n\n");  /* Show cursor again */
    return Particles;
}


/* -------------------------
 * Functions of phase 3: Visualization
 ------------------------- */

RGB *create_image(int width, int height) {
    RGB *image;
    int pixel_count, i;
    pixel_count = (size_t)width * (size_t)height;
    image = malloc(pixel_count * sizeof(RGB));


    if (image == NULL) {
        return NULL;
    }
    /* Store pixels contiguously in row-major order and initialize the canvas
    to black before trajectories and nuclei are drawn.
     */
    for (i = 0; i < pixel_count; i++) {
        image[i].r = 0;
        image[i].g = 0;
        image[i].b = 0;
    }

    return image;
}


/*Set one image pixel if its coordinates are in bounds. */
Validate color_pixel(int image_width, int image_height, int iy, int ix, RGB *image, RGB color) {
    if (ix < 0 || ix >= image_width || iy < 0 || iy >= image_height) {
        return INVALID; /* Nothing to draw inside the image */
    }
    int position = iy * image_width + ix;
    image[position].r = color.r;
    image[position].g = color.g;
    image[position].b = color.b;
    return VALID;
}

/* Rasterize the segment between two physical positions using evenly spaced
samples (a digital differential analyzer). Pixel samples outside the image
are ignored; only in-bounds pixels are written.
 */
Validate draw_trajectories(Vector2 prev_pos, Vector2 pos, RGB color, double physical_width, int image_width, double physical_height, int image_height, RGB *image) {
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
        if (pp.x >= 0 && pp.x < image_width && pp.y >= 0 && pp.y < image_height) {
            color_pixel(image_width, image_height, pp.y, pp.x, image, color);
        }
        return VALID;
    }

    double x_step = (double)dx / steps;
    double y_step = (double)dy / steps;

    for (int i = 0; i <= steps; i++) {
        ix = (int)round(x);
        iy = (int)round(y);

        if (ix >= 0 && ix < image_width && iy >= 0 && iy < image_height) {
            if (color_pixel(image_width, image_height, iy, ix, image, color) == INVALID) {
                return INVALID;
            }
        }

        x += x_step;
        y += y_step;
    }

    return VALID;

}


void draw_nuclei(Nucleus *nuclei,int nucleus_count,SimulationParameters SP,RGB *image){
    for (int i = 0; i < nucleus_count; i++){
        draw_single_nucleus(nuclei[i],SP,image);
    }
    
    return;
}


void draw_single_nucleus(Nucleus nuclei, SimulationParameters SP,RGB *image){
    //set r 1-max pixels now 9 based on ratios
    int r = (int)round(2.0 + 0.8 * sqrt((double)nuclei.Z));

    int sim_x = physical_to_pixel_x(nuclei.pos.x,SP.physical_width, SP.image_width);
    int sim_y = physical_to_pixel_y(nuclei.pos.y,SP.physical_height, SP.image_height);
    RGB color = {255,255,255};
    
    for (int i = sim_x-r; i <= sim_x+r; i++){
        int dx = abs(sim_x-i);
        for (int j = sim_y-r; j <= sim_y+r; j++){
            int dy = abs(sim_y-j);
            if (dx*dx + dy*dy <= r*r && i >= 0 && i < SP.image_width && j >= 0 && j < SP.image_height) {
                color_pixel(SP.image_width, SP.image_height, j, i, image, color);
            }
        }
    }
    
    return;
}


/* Display the RGB buffer in an SDL 1.2 window and keep it open until the user
closes the window. SDL surfaces created here are released on exit/errors.
 */



Validate draw(const RGB *image, int width, int height) { // Takes image and converts SDL compatible and draws
    SDL_Event ev; //Create the exit
    if (SDL_Init(SDL_INIT_VIDEO) < 0) { //Initialises SDL
        fprintf(stderr, "SDL init error: %s\n", SDL_GetError());
        return INVALID;
    }
    SDL_Surface *screen = SDL_SetVideoMode(width, height, 32, SDL_SWSURFACE); // screen will be shown
    if (screen == NULL) {
        fprintf(stderr, "SDL screen error: %s\n", SDL_GetError());
        SDL_Quit();
        return INVALID;
    }
    SDL_Surface *surface = SDL_CreateRGBSurface(SDL_SWSURFACE, width, height, 32, 0, 0, 0, 0); //SDL Compatible data storage
    if (surface == NULL) {
        fprintf(stderr, "SDL error: %s\n", SDL_GetError());
        SDL_Quit();
        return INVALID;
    }
    if (SDL_MUSTLOCK(surface) && SDL_LockSurface(surface) < 0) {
        fprintf(stderr, "SDL lock error: %s\n", SDL_GetError());
        SDL_FreeSurface(surface);
        SDL_Quit();
        return INVALID;
    }
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            Uint32 color = SDL_MapRGB(surface->format,image[width * i + j].r,image[width * i + j].g,image[width * i + j].b);
            Uint32 *pixel = (Uint32 *)((Uint8 *)surface->pixels + i * surface->pitch) + j;
            *pixel = color;
        }
    }
    if (SDL_MUSTLOCK(surface)) {
        SDL_UnlockSurface(surface);
    }
    if (SDL_BlitSurface(surface, NULL, screen, NULL) < 0) {
        fprintf(stderr, "SDL blit error: %s\n", SDL_GetError());
        SDL_FreeSurface(surface);
        SDL_Quit();
        return INVALID;
    }
    if (SDL_Flip(screen) < 0) {
        fprintf(stderr, "SDL flip error: %s\n", SDL_GetError());
        SDL_FreeSurface(surface);
        SDL_Quit();
        return INVALID;
    }
    while (SDL_WaitEvent(&ev) && ev.type != SDL_QUIT) {}
    SDL_FreeSurface(surface);
    SDL_Quit();
    return VALID;
}


Validate save_ppm(const char *filename, const RGB *image, int width, int height) {
    FILE *ppm=fopen(filename, "w");
    if (ppm == NULL) {
        fprintf(stderr,"Unable to open file for writing.\n");
        return INVALID;
    }
    fprintf(ppm,"P3\n %d %d\n 255\n", width, height);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            fprintf(ppm, "%d %d %d ", image[i*width+j].r,image[i*width+j].g,image[i*width+j].b);
        }
        fprintf(ppm, "\n");
    }
    fclose(ppm);
    return VALID;
}


/* -------------------------
 * Main: Program modes
 ------------------------- */

/* Run the interactive simulation path: collect settings, load nuclei, create
and simulate particles, save/display the image, then release owned buffers.
 */
int normal_mode(void) {
    int nucleus_count = 0; /*  Set by read_nuclei after parsing the input file. */

    Nucleus *nuclei = NULL;
    AlphaParticle *Particles = NULL; /* Owned particle array, released before return. */
    SimulationParameters SIM_PAR;

    if (read_parameters(&SIM_PAR) == INVALID) {
        fprintf(stderr, "Invalid input. Returning to main menu.\n");
        return 1;
    }
    if (validate_parameters(SIM_PAR) == VALID){
        printf("Simulation parameters validated.\n");
    }
    else{
        fprintf(stderr,"There is one or more invalid parameters.\n");
        return 1;
    }

    if (validate_file(SIM_PAR) == INVALID) {
        return 1;
    }

    RGB *image = create_image(SIM_PAR.image_width, SIM_PAR.image_height);
    if (image == NULL) {
        fprintf(stderr,"Memory allocation for image failed\n");
        return 1;
    }
    nuclei = read_nuclei(SIM_PAR, &nucleus_count);
    if (nuclei == NULL && nucleus_count == 0) {
        fprintf(stderr,"something went wrong with the file.\n");
        free(image);
        return 1;
    }

    if (run_symulation(Particles, SIM_PAR, image, nuclei, nucleus_count)) {
        return 0;
    }
    return 1;
}


int demo_mode(void) {
    SimulationParameters SIM_PAR = read_parameters_demo();
    int nucleus_count = 7;


    AlphaParticle *Particles = NULL; /* Owned particle array, released before return. */

    RGB *image = create_image(SIM_PAR.image_width, SIM_PAR.image_height);
    if (image == NULL) {
        fprintf(stderr,"Memory allocation for image failed\n");
        return 1;
    }

    Nucleus *nuclei = malloc(nucleus_count * sizeof(Nucleus));

    if (nuclei == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(image);
        return 1;
    }

    Nucleus initial_nuclei[] = {
        {{5.3e-13, 1.2e-13}, 35},
        {{7.2e-13, 1.7e-13}, 60},
        {{8.6e-13, 3.1e-13}, 79},
        {{9.1e-13, 5.0e-13}, 45},
        {{8.6e-13, 6.9e-13}, 55},
        {{7.2e-13, 8.3e-13}, 79},
        {{5.3e-13, 8.8e-13}, 85}
    };

    memcpy(nuclei, initial_nuclei, sizeof(initial_nuclei));

    if (run_symulation(Particles, SIM_PAR, image, nuclei, nucleus_count)) {
        return 0;
    }
    return 1;
}


int read_mode(void) {
    SimulationParameters SIM_PAR;
    int nucleus_count = 0; /* Set by read_nuclei after parsing the input file. */

    Nucleus *nuclei = NULL;
    AlphaParticle *Particles = NULL; /* Owned particle array, released before return. */

    if (read_parameters_from_file(&SIM_PAR) == INVALID) {
        fprintf(stderr, "Invalid input. Returning to main menu.\n");
        return 1;
    }
    if (validate_parameters(SIM_PAR) == VALID){
        printf("Simulation parameters validated.\n");
    }
    else{
        fprintf(stderr,"There is one or more invalid parameters.\n");
        return 1;
    }

    if (validate_file(SIM_PAR) == INVALID) {
        return 1;
    }

    RGB *image = create_image(SIM_PAR.image_width, SIM_PAR.image_height);
    if (image == NULL) {
        fprintf(stderr,"Memory allocation for image failed\n");
        return 1;
    }
    nuclei = read_nuclei(SIM_PAR, &nucleus_count);
    if (nuclei == NULL && nucleus_count == 0) {
        fprintf(stderr,"something went wrong with the file.\n");
        free(image);
        return 1;
    }

    if (run_symulation(Particles, SIM_PAR, image, nuclei, nucleus_count)) {
        return 0;
    }
    return 1;
}


Validate run_symulation(AlphaParticle *Particles, SimulationParameters SIM_PAR, RGB *image, Nucleus *nuclei, int nucleus_count) {
    /* Initialize and integrate the particle array, adding paths to image. */
    Particles = initialize_particles(SIM_PAR);
    if (Particles == NULL) {
        fprintf(stderr,"something went wrong with with initialising the particles\n");
        free(image);
        free(nuclei);
        return INVALID;
    }
    if (simulate_particles(nuclei, Particles, nucleus_count, SIM_PAR, image) == NULL) {
        fprintf(stderr,"simulate_particles failed\n");
        free(image);
        free(nuclei);
        free(Particles);
        return INVALID;
    }

    draw_nuclei(nuclei, nucleus_count, SIM_PAR, image);

    if (draw(image, SIM_PAR.image_width, SIM_PAR.image_height) == INVALID) {
        free(image);
        free(nuclei);
        free(Particles);
        return INVALID;
    }

    printf("Do you want to save the image as a PPM file? (y/n): ");
    char save_choice;
    if (scanf(" %c", &save_choice) == 1 && save_choice == 'y') {
        if (save_ppm(SIM_PAR.output_filename, image, SIM_PAR.image_width, SIM_PAR.image_height) == INVALID) {
            fprintf(stderr, "save_ppm failed\n");
            free(image);
            free(nuclei);
            free(Particles);
            return INVALID;
        }
        printf("Image saved as %s\n", SIM_PAR.output_filename);
    }
    free(image);
    free(Particles);
    free(nuclei);
    return VALID;
}


/* -------------------------
* Main
------------------------- */


    int main(void) {
        srand(time(NULL)); /* for seed generation */
        int choice, True = 0;
        /* Keep presenting the menu until the user selects the exit option. */
        while (True == 0) {
            printf("\nSelect mode:\n- 1:Demo mode\n- 2:Read simulation parameters from file\n- 3:Enter simulation parameters by hand\n- 4:HELP\n- 9:EXIT\n");
            if (scanf("%d", &choice) != 1) {
                fprintf(stderr, "Invalid input.\n");
                while (getchar() != '\n'); /* Clear the input buffer */
                continue; /* Resets the loop to present the menu again */
            }
            switch (choice) {
                case 1:
                    demo_mode();
                    break;
                case 2:
                    read_mode();
                    break;
                case 3:
                    normal_mode();
                    break;
                case 4:
                    printf("\n========== RUTHERFORD SCATTERING SIMULATION ==========\n\n");

                    printf("This program simulates alpha particles travelling through\n");
                    printf("a field of positively charged atomic nuclei.\n");
                    printf("Their trajectories are calculated using Coulomb's law\n");
                    printf("and the fourth-order Runge-Kutta (RK4) method.\n\n");

                    printf("MODES:\n");
                    printf("1 - Demo: Run a simulation with predefined parameters and nuclei.\n");
                    printf("2 - File: Load simulation parameters and nuclei from TXT files.\n");
                    printf("3 - Manual: Enter parameters manually and load nuclei from a TXT file.\n");
                    printf("9 - Exit the program.\n\n");

                    printf("SIMULATION PARAMETERS:\n");
                    printf("- Image width and height: Output resolution in pixels.\n");
                    printf("- Physical width and height: Simulation area in metres.\n");
                    printf("- Particle count: Number of alpha particles.\n");
                    printf("- Starting X: Initial horizontal position in metres.\n");
                    printf("- Velocity min/max: Initial velocity range in m/s.\n");
                    printf("- dt: Time interval between numerical integration steps.\n");
                    printf("- Max steps: Number of integration steps per particle.\n\n");

                    printf("NUCLEI FILE FORMAT (.txt):\n");
                    printf("Each line contains: X Y Z\n");
                    printf("Example: 5.3e-13 1.2e-13 79\n");
                    printf("Z is optional (default: 79).\n");
                    printf("Allowed atomic numbers: 1-92.\n");
                    printf("Nuclei must be inside the simulation area,\n");
                    printf("in the right half of the screen.\n\n");

                    printf("PARAMETER FILE FORMAT (.txt):\n");
                    printf("Ten whitespace-separated values in this order:\n");
                    printf("image_width image_height physical_width physical_height\n");
                    printf("particle_count starting_x velocity_min velocity_max dt max_steps\n\n");

                    printf("OUTPUT:\n");
                    printf("The trajectories are displayed in an SDL window.\n");
                    printf("Close the window to return to the terminal.\n");
                    printf("You can then save the image in PPM format.\n\n");

                    printf("=====================================================\n\n");
                    break;
                case 9:
                    True = 1; /* Exit the loop and terminate the program */
                    break;
                default:
                    printf("Invalid choice.\n");
                    break;
            }
        }
    }
