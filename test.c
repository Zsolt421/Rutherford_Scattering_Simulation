#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef enum{
    INVALID,
    VALID,
} Validate;

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

Validate validate_nuclei(SimulationParameters SIMPAR, int read_count, int count, double x, double y, int Z){
    if (x > SIMPAR.physical_width || y > SIMPAR.physical_height){
        fprintf(stderr, "Nuclei position invalid in %s line %d", SIMPAR.nuclei_filename, count+1);
        return INVALID;
    }

    if (read_count == 2){
        return VALID;
    }

    if (Z <= 0 || Z > 92){ //Atomic number check, exists but not larger than Uranium
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

    file = fopen(SIMPAR.nuclei_filename, "r");
    if (file == NULL) {
        *out_count = 0;
        return NULL;
    }

    while (fgets(line, sizeof(line), file) != NULL) {

        read_count = sscanf(line, "%lf %lf %d", &x, &y, &Z);

        if (read_count >= 2) {

            if (validate_nuclei(SIMPAR,read_count, count,x,y,Z) == INVALID){
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

            count++;
        }
        
        // new part !!
        else {
            fprintf(stderr, "Not enough parameters in %s, line %d", SIMPAR.nuclei_filename, count+1);
            return NULL;
        }
        // new part !!
    }

    fclose(file);
    *out_count = count;
    return p;
}


int main() {
    int nucleus_count = 0; /* will be scanned from file, 0 is only placeholder */
    Nucleus *nuclei = NULL;
    SimulationParameters SIM_PAR = {100, 1000, 1000, 1000,1000,1000,1000,1000,1000, "data.txt", "alma.txt"};

    nuclei = read_nuclei(SIM_PAR, &nucleus_count); 
    if (nuclei == NULL) return 1;

    for (int i =0; i < nucleus_count; i++){
        printf("Nuclei %d: %lf %lf %d\n", i+1, nuclei[i].x,nuclei[i].y,nuclei[i].Z);
    }
    free(nuclei);
}