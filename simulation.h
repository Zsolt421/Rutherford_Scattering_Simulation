#ifndef SIMULATION_H
#define SIMULATION_H

/* -------------------------
 * Definition of structures
 * ------------------------- */

typedef enum {
    INVALID,
    VALID
} Validate;

typedef struct {
    double x;
    double y;
} Vector2;

typedef struct {
    Vector2 pos;
    int Z;
} Nucleus;

typedef struct {
    unsigned char r;
    unsigned char g;
    unsigned char b;
} RGB;

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

/* -------------------------
 * Definition of constants
 * ------------------------- */

#define PI 3.14159265358979
#define EPSILON_0 8.8541878188e-12
#define COULOMB_CONSTANT (1.0 / (4.0 * PI * EPSILON_0))
#define ELEMENTARY_CHARGE 1.602176634e-19
#define ALPHA_CHARGE (2.0 * ELEMENTARY_CHARGE)
#define ALPHA_MASS 6.644657e-27

/* -------------------------
 * General helper functions
 * ------------------------- */

Vector2 vec_add(Vector2 a, Vector2 b);
Vector2 vec_scale(Vector2 v, double factor);
int physical_to_pixel_x(double x, double physical_width, int image_width);
int physical_to_pixel_y(double y, double physical_height, int image_height);
void show_progress(long long completed, long long total);

/* -------------------------
 * Phase 2: Input and validation
 * ------------------------- */

char *txt_adder(char *filename);
char *ppm_adder(char *filename);
SimulationParameters read_parameters_demo(void);
Validate read_parameters(SimulationParameters *SIMPAR);
Validate read_parameters_from_file(SimulationParameters *SIMPAR);
Validate validate_parameters(SimulationParameters SIMPAR);
Validate validate_file(SimulationParameters SP);
Validate validate_nuclei(SimulationParameters SIMPAR, int read_count, int count, double x, double y, int Z);
Nucleus *read_nuclei(SimulationParameters SIMPAR, int *out_count);

/* -------------------------
 * Phase 1: Physics and simulation
 * ------------------------- */

Vector2 calculate_coulomb_force(Vector2 position, Nucleus nucleus);
Vector2 calculate_total_force(Vector2 position, const Nucleus nuclei[], int nucleus_count);
Validate numerical_calculations(AlphaParticle *alpha, const Nucleus nuclei[], int nucleus_count, SimulationParameters SIMPAR, RGB *image, int particle_index);
AlphaParticle *initialize_particles(SimulationParameters SIMPAR);
AlphaParticle *simulate_particles(const Nucleus nuclei[], AlphaParticle *Particles, int nucleus_count, SimulationParameters SIMPAR, RGB *image);

/* -------------------------
 * Phase 3: Visualization
 * ------------------------- */

RGB *create_image(int width, int height);
Validate color_pixel(int image_width, int image_height, int iy, int ix, RGB *image, RGB color);
Validate draw_trajectories(Vector2 prev_pos, Vector2 pos, RGB color, double physical_width, int image_width, double physical_height, int image_height, RGB *image);
void draw_nuclei(Nucleus *nuclei, int nucleus_count, SimulationParameters SP, RGB *image);
void draw_single_nucleus(Nucleus nucleus, SimulationParameters SP, RGB *image);
Validate draw(const RGB *image, int width, int height);
Validate save_ppm(const char *filename, const RGB *image, int width, int height);

/* -------------------------
 * Main: Program modes
 * ------------------------- */

int normal_mode(void);
int demo_mode(void);
int read_mode(void);
Validate run_symulation(AlphaParticle *Particles, SimulationParameters SIM_PAR, RGB *image, Nucleus *nuclei, int nucleus_count);

#endif /* SIMULATION_H */
