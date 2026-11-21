#ifndef SIMULATION_H
#define SIMULATION_H


/** -------------------------
* Definition of structures
 ------------------------- */

typedef enum{
    INVALID,
    VALID,
} Validate;
/** Can be used to store position(x, y), velocity(vx, vy), force(Fx, Fy), acceleration(ax,ay). */
typedef struct
{
    double x;
    double y;
} Vector2;
/** Position and atomic number (positive charge, defaults to 79 if missing). */
typedef struct {
    Vector2 pos;
    int Z;
} Nucleus;
/** Color definition for unique colors. */
typedef struct {
    unsigned char r;
    unsigned char g;
    unsigned char b;
} RGB;
/** Stores data for each alpha particle. */
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
    double initial_velocity_min;
    double initial_velocity_max;
    double dt;
    int max_steps;
    char nuclei_filename[256];
    char output_filename[256];
} SimulationParameters;

typedef struct {
    int x;
    int y;
} Pixel2;


/** -------------------------
* Definition of constants
 ------------------------- */
#define PI 3.14159265358979
#define EPSILON_0 8.8541878188e-12
#define COULOMB_CONSTANT (1.0 / (4.0 * PI * EPSILON_0))

#define ELEMENTARY_CHARGE 1.602176634e-19
#define ALPHA_CHARGE (2.0 * ELEMENTARY_CHARGE)
#define ALPHA_MASS 6.644657e-27


/** -------------------------
* Definition of prototypes
 ------------------------- */


/** Helper functions */
Vector2 vec_add(Vector2 a, Vector2 b);
Vector2 vec_scale(Vector2 v, double factor);
int physical_to_pixel_x(double x, double physical_width, int image_width);
int physical_to_pixel_y(double y, double physical_height, int image_height);

/** Phase 2: Input and validation */
SimulationParameters read_parameters(void);
Validate validate_parameters(SimulationParameters SIMPAR);
Validate validate_file(SimulationParameters SP);
Validate validate_nuclei(SimulationParameters SIMPAR, int read_count, int count, double x, double y, int Z);
Nucleus *read_nuclei(SimulationParameters SIMPAR, int *out_count);

/** Phase 1: Physics and simulation */
Vector2 calculate_coulomb_force(Vector2 position, Nucleus nucleus);
Vector2 calculate_total_force(Vector2 position, const Nucleus nuclei[], int nucleus_count);
Validate numerical_calculations(AlphaParticle *alpha, const Nucleus nuclei[], int nucleus_count, SimulationParameters SIMPAR, RGB *image);
AlphaParticle *initialize_particles(SimulationParameters SIMPAR);
AlphaParticle *simulate_particles(const Nucleus nuclei[], AlphaParticle *Particles, int nucleus_count, SimulationParameters SIMPAR, RGB *image);

/** Phase 3: Visualization */
RGB *create_image(int width, int height);
Validate color_pixel(int image_width, int image_height, int iy, int ix, RGB *image, RGB color);
Validate draw_trajectories(Vector2 prev_pos, Vector2 pos, RGB color, double physical_width, int image_width, double physical_height, int image_height, RGB *image);
void draw_nuclei(Nucleus *nuclei, int nucleus_count, int r, SimulationParameters SP, RGB *image);
void draw_single_nucleus(Nucleus nucleus, int r, SimulationParameters SP, RGB *image);
Validate draw(const RGB *image, int width, int height);
Validate save_ppm(const char *filename, const RGB *image, int width, int height);


#endif

#ifndef RUTHERFORD_SCATTERING_SIMULATION_SIMULATION_H
#define RUTHERFORD_SCATTERING_SIMULATION_SIMULATION_H

/** RUTHERFORD_SCATTERING_SIMULATION_SIMULATION_H */
#endif
