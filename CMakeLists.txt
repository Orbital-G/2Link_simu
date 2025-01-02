cmake_minimum_required(VERSION 3.10)

project(2LinkSimu)
set(CMAKE_CXX_STANDARD 17)

add_executable(simu_with_drawing
    src/simu_with_drawing.C
)

target_include_directories(simu_with_drawing PRIVATE
    ${PROJECT_SOURCE_DIR}/include

    /usr/local/Cellar/eigen/3.4.0_1/include/eigen3
)

# SFML ライブラリをリンク
find_package(SFML 2.5 COMPONENTS graphics window system REQUIRED)
target_link_libraries(simu_with_drawing PRIVATE sfml-graphics sfml-window sfml-system)
