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

void draw(const RGB *image, int width, int height) { // Takes image and converts SDL compatible and draws
    SDL_Event ev; //Create the exit
    if (SDL_Init(SDL_INIT_VIDEO) < 0) { //Initialises SDL
        fprintf(stderr, "SDL init error: %s\n", SDL_GetError());
        return;
    }
    SDL_Surface *screen = SDL_SetVideoMode(width, height, 32, SDL_SWSURFACE); // screen will be shown
    if (screen == NULL) {
        fprintf(stderr, "SDL screen error: %s\n", SDL_GetError());
        SDL_Quit();
        return;
    }
    SDL_Surface *surface = SDL_CreateRGBSurface(SDL_SWSURFACE, width, height, 32, 0, 0, 0, 0); //SDL Compatible data storage
    if (surface == NULL) {
        fprintf(stderr, "SDL error: %s\n", SDL_GetError());
        SDL_Quit();
        return;
    }
    if (SDL_MUSTLOCK(surface) && SDL_LockSurface(surface) < 0) {
        fprintf(stderr, "SDL lock error: %s\n", SDL_GetError());
        SDL_FreeSurface(surface);
        SDL_Quit();
        return;
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
        return;
    }
    if (SDL_Flip(screen) < 0) {
        fprintf(stderr, "SDL flip error: %s\n", SDL_GetError());
        SDL_FreeSurface(surface);
        SDL_Quit();
        return;
    }
    while (SDL_WaitEvent(&ev) && ev.type != SDL_QUIT) {}
    SDL_FreeSurface(surface);
    SDL_Quit();
}
