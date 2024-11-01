cmake_minimum_required(VERSION 3.30 FATAL_ERROR)

project(gasdyn LANGUAGES Fortran)

set(gasdyn-libdir ${CMAKE_BINARY_DIR}/lib)
set(gasdyn-bindir ${CMAKE_BINARY_DIR}/bin)

add_executable(${PROJECT_NAME})

set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES Fortran_STANDARD 2018
               Fortran_STANDARD_REQUIRED on
               Fortran_MODULE_DIRECTORY ${gasdyn-bindir}
               ARCHIVE_OUTPUT_DIRECTORY ${gasdyn-libdir}
               LIBRARY_OUTPUT_DIRECTORY ${gasdyn-libdir}
               RUNTIME_OUTPUT_DIRECTORY ${gasdyn-bindir}
)

add_subdirectory(source)

add_subdirectory(input)

add_subdirectory(plots)
