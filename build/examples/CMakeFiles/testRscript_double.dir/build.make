# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.19

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
CMAKE_COMMAND = /lambda_stor/homes/xyu/cmake/cmake-3.19.1-Linux-x86_64/bin/cmake

# The command to remove a file.
RM = /lambda_stor/homes/xyu/cmake/cmake-3.19.1-Linux-x86_64/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /homes/xyu/cuZC/Z-checker

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /homes/xyu/cuZC/Z-checker/build

# Include any dependencies generated for this target.
include examples/CMakeFiles/testRscript_double.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/testRscript_double.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/testRscript_double.dir/flags.make

examples/CMakeFiles/testRscript_double.dir/testRscript_double.c.o: examples/CMakeFiles/testRscript_double.dir/flags.make
examples/CMakeFiles/testRscript_double.dir/testRscript_double.c.o: ../examples/testRscript_double.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/homes/xyu/cuZC/Z-checker/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object examples/CMakeFiles/testRscript_double.dir/testRscript_double.c.o"
	cd /homes/xyu/cuZC/Z-checker/build/examples && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/testRscript_double.dir/testRscript_double.c.o -c /homes/xyu/cuZC/Z-checker/examples/testRscript_double.c

examples/CMakeFiles/testRscript_double.dir/testRscript_double.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/testRscript_double.dir/testRscript_double.c.i"
	cd /homes/xyu/cuZC/Z-checker/build/examples && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /homes/xyu/cuZC/Z-checker/examples/testRscript_double.c > CMakeFiles/testRscript_double.dir/testRscript_double.c.i

examples/CMakeFiles/testRscript_double.dir/testRscript_double.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/testRscript_double.dir/testRscript_double.c.s"
	cd /homes/xyu/cuZC/Z-checker/build/examples && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /homes/xyu/cuZC/Z-checker/examples/testRscript_double.c -o CMakeFiles/testRscript_double.dir/testRscript_double.c.s

# Object files for target testRscript_double
testRscript_double_OBJECTS = \
"CMakeFiles/testRscript_double.dir/testRscript_double.c.o"

# External object files for target testRscript_double
testRscript_double_EXTERNAL_OBJECTS =

bin/testRscript_double: examples/CMakeFiles/testRscript_double.dir/testRscript_double.c.o
bin/testRscript_double: examples/CMakeFiles/testRscript_double.dir/build.make
bin/testRscript_double: lib/libzc.so
bin/testRscript_double: examples/CMakeFiles/testRscript_double.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/homes/xyu/cuZC/Z-checker/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable ../bin/testRscript_double"
	cd /homes/xyu/cuZC/Z-checker/build/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testRscript_double.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/testRscript_double.dir/build: bin/testRscript_double

.PHONY : examples/CMakeFiles/testRscript_double.dir/build

examples/CMakeFiles/testRscript_double.dir/clean:
	cd /homes/xyu/cuZC/Z-checker/build/examples && $(CMAKE_COMMAND) -P CMakeFiles/testRscript_double.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/testRscript_double.dir/clean

examples/CMakeFiles/testRscript_double.dir/depend:
	cd /homes/xyu/cuZC/Z-checker/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /homes/xyu/cuZC/Z-checker /homes/xyu/cuZC/Z-checker/examples /homes/xyu/cuZC/Z-checker/build /homes/xyu/cuZC/Z-checker/build/examples /homes/xyu/cuZC/Z-checker/build/examples/CMakeFiles/testRscript_double.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/testRscript_double.dir/depend

