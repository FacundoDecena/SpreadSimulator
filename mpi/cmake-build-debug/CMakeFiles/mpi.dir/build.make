# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/facundo/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/201.7846.88/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/facundo/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/201.7846.88/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/facundo/Datos/LoQueUso/UNSL/2020/Primer-Cuatrimestre/Paralelos/PracticoFinal/COVID-19/mpi

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/facundo/Datos/LoQueUso/UNSL/2020/Primer-Cuatrimestre/Paralelos/PracticoFinal/COVID-19/mpi/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/mpi.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mpi.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mpi.dir/flags.make

CMakeFiles/mpi.dir/main.c.o: CMakeFiles/mpi.dir/flags.make
CMakeFiles/mpi.dir/main.c.o: ../main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/facundo/Datos/LoQueUso/UNSL/2020/Primer-Cuatrimestre/Paralelos/PracticoFinal/COVID-19/mpi/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/mpi.dir/main.c.o"
	/usr/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/mpi.dir/main.c.o   -c /media/facundo/Datos/LoQueUso/UNSL/2020/Primer-Cuatrimestre/Paralelos/PracticoFinal/COVID-19/mpi/main.c

CMakeFiles/mpi.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mpi.dir/main.c.i"
	/usr/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /media/facundo/Datos/LoQueUso/UNSL/2020/Primer-Cuatrimestre/Paralelos/PracticoFinal/COVID-19/mpi/main.c > CMakeFiles/mpi.dir/main.c.i

CMakeFiles/mpi.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mpi.dir/main.c.s"
	/usr/bin/mpicc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /media/facundo/Datos/LoQueUso/UNSL/2020/Primer-Cuatrimestre/Paralelos/PracticoFinal/COVID-19/mpi/main.c -o CMakeFiles/mpi.dir/main.c.s

# Object files for target mpi
mpi_OBJECTS = \
"CMakeFiles/mpi.dir/main.c.o"

# External object files for target mpi
mpi_EXTERNAL_OBJECTS =

mpi: CMakeFiles/mpi.dir/main.c.o
mpi: CMakeFiles/mpi.dir/build.make
mpi: CMakeFiles/mpi.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/facundo/Datos/LoQueUso/UNSL/2020/Primer-Cuatrimestre/Paralelos/PracticoFinal/COVID-19/mpi/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable mpi"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mpi.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mpi.dir/build: mpi

.PHONY : CMakeFiles/mpi.dir/build

CMakeFiles/mpi.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mpi.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mpi.dir/clean

CMakeFiles/mpi.dir/depend:
	cd /media/facundo/Datos/LoQueUso/UNSL/2020/Primer-Cuatrimestre/Paralelos/PracticoFinal/COVID-19/mpi/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/facundo/Datos/LoQueUso/UNSL/2020/Primer-Cuatrimestre/Paralelos/PracticoFinal/COVID-19/mpi /media/facundo/Datos/LoQueUso/UNSL/2020/Primer-Cuatrimestre/Paralelos/PracticoFinal/COVID-19/mpi /media/facundo/Datos/LoQueUso/UNSL/2020/Primer-Cuatrimestre/Paralelos/PracticoFinal/COVID-19/mpi/cmake-build-debug /media/facundo/Datos/LoQueUso/UNSL/2020/Primer-Cuatrimestre/Paralelos/PracticoFinal/COVID-19/mpi/cmake-build-debug /media/facundo/Datos/LoQueUso/UNSL/2020/Primer-Cuatrimestre/Paralelos/PracticoFinal/COVID-19/mpi/cmake-build-debug/CMakeFiles/mpi.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mpi.dir/depend
