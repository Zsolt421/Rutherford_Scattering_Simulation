cmake_minimum_required(VERSION 3.16)

project(Rutherford_Scattering_Simulation C)

add_executable(Rutherford main.c)
add_executable(Szoszmoty Szoszmoty.c)
# SDL 1.2 helye a sajat gepemen
set(SDL_ROOT
        "C:/Users/István/OneDrive - movarinfo/Asztal/Studies/Programozás/BHW/mingw-w64-ucrt-x86_64-SDL-1.2.15+r419+gef3a6c05-2-any.pkg/ucrt64"
)

# Az SDL fejlecfajljainak eleresi utja
target_include_directories(Rutherford PRIVATE
        "${SDL_ROOT}/include"
)

target_link_libraries(Rutherford PRIVATE
        mingw32
        "${SDL_ROOT}/lib/libSDLmain.a"
        "${SDL_ROOT}/lib/libSDL.dll.a"
)
# SDL.dll masolasa a program melle
add_custom_command(TARGET Rutherford POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${SDL_ROOT}/bin/SDL.dll"
        "$<TARGET_FILE_DIR:Rutherford>"
)