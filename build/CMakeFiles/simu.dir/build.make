# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.31

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/kcbcp/Labo_document/2Link_simu

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/kcbcp/Labo_document/2Link_simu/build

# Include any dependencies generated for this target.
include CMakeFiles/simu.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/simu.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/simu.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/simu.dir/flags.make

CMakeFiles/simu.dir/codegen:
.PHONY : CMakeFiles/simu.dir/codegen

CMakeFiles/simu.dir/src/simu_with_drawing.C.o: CMakeFiles/simu.dir/flags.make
CMakeFiles/simu.dir/src/simu_with_drawing.C.o: /Users/kcbcp/Labo_document/2Link_simu/src/simu_with_drawing.C
CMakeFiles/simu.dir/src/simu_with_drawing.C.o: CMakeFiles/simu.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/kcbcp/Labo_document/2Link_simu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/simu.dir/src/simu_with_drawing.C.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/simu.dir/src/simu_with_drawing.C.o -MF CMakeFiles/simu.dir/src/simu_with_drawing.C.o.d -o CMakeFiles/simu.dir/src/simu_with_drawing.C.o -c /Users/kcbcp/Labo_document/2Link_simu/src/simu_with_drawing.C

CMakeFiles/simu.dir/src/simu_with_drawing.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/simu.dir/src/simu_with_drawing.C.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/kcbcp/Labo_document/2Link_simu/src/simu_with_drawing.C > CMakeFiles/simu.dir/src/simu_with_drawing.C.i

CMakeFiles/simu.dir/src/simu_with_drawing.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/simu.dir/src/simu_with_drawing.C.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/kcbcp/Labo_document/2Link_simu/src/simu_with_drawing.C -o CMakeFiles/simu.dir/src/simu_with_drawing.C.s

CMakeFiles/simu.dir/include/linearAlgebra/vector/vector.C.o: CMakeFiles/simu.dir/flags.make
CMakeFiles/simu.dir/include/linearAlgebra/vector/vector.C.o: /Users/kcbcp/Labo_document/2Link_simu/include/linearAlgebra/vector/vector.C
CMakeFiles/simu.dir/include/linearAlgebra/vector/vector.C.o: CMakeFiles/simu.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/kcbcp/Labo_document/2Link_simu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/simu.dir/include/linearAlgebra/vector/vector.C.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/simu.dir/include/linearAlgebra/vector/vector.C.o -MF CMakeFiles/simu.dir/include/linearAlgebra/vector/vector.C.o.d -o CMakeFiles/simu.dir/include/linearAlgebra/vector/vector.C.o -c /Users/kcbcp/Labo_document/2Link_simu/include/linearAlgebra/vector/vector.C

CMakeFiles/simu.dir/include/linearAlgebra/vector/vector.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/simu.dir/include/linearAlgebra/vector/vector.C.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/kcbcp/Labo_document/2Link_simu/include/linearAlgebra/vector/vector.C > CMakeFiles/simu.dir/include/linearAlgebra/vector/vector.C.i

CMakeFiles/simu.dir/include/linearAlgebra/vector/vector.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/simu.dir/include/linearAlgebra/vector/vector.C.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/kcbcp/Labo_document/2Link_simu/include/linearAlgebra/vector/vector.C -o CMakeFiles/simu.dir/include/linearAlgebra/vector/vector.C.s

# Object files for target simu
simu_OBJECTS = \
"CMakeFiles/simu.dir/src/simu_with_drawing.C.o" \
"CMakeFiles/simu.dir/include/linearAlgebra/vector/vector.C.o"

# External object files for target simu
simu_EXTERNAL_OBJECTS =

simu: CMakeFiles/simu.dir/src/simu_with_drawing.C.o
simu: CMakeFiles/simu.dir/include/linearAlgebra/vector/vector.C.o
simu: CMakeFiles/simu.dir/build.make
simu: /usr/local/lib/libsfml-graphics.2.6.2.dylib
simu: /usr/local/lib/libsfml-window.2.6.2.dylib
simu: /usr/local/lib/libsfml-system.2.6.2.dylib
simu: CMakeFiles/simu.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/kcbcp/Labo_document/2Link_simu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable simu"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/simu.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/simu.dir/build: simu
.PHONY : CMakeFiles/simu.dir/build

CMakeFiles/simu.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/simu.dir/cmake_clean.cmake
.PHONY : CMakeFiles/simu.dir/clean

CMakeFiles/simu.dir/depend:
	cd /Users/kcbcp/Labo_document/2Link_simu/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/kcbcp/Labo_document/2Link_simu /Users/kcbcp/Labo_document/2Link_simu /Users/kcbcp/Labo_document/2Link_simu/build /Users/kcbcp/Labo_document/2Link_simu/build /Users/kcbcp/Labo_document/2Link_simu/build/CMakeFiles/simu.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/simu.dir/depend
