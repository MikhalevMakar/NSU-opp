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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /tmp/tmp.v25c7svKbI

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /tmp/tmp.v25c7svKbI/cmake-build-debug-remote-host-1

# Include any dependencies generated for this target.
include CMakeFiles/Version1.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Version1.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Version1.dir/flags.make

CMakeFiles/Version1.dir/main.cpp.o: CMakeFiles/Version1.dir/flags.make
CMakeFiles/Version1.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/tmp/tmp.v25c7svKbI/cmake-build-debug-remote-host-1/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Version1.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Version1.dir/main.cpp.o -c /tmp/tmp.v25c7svKbI/main.cpp

CMakeFiles/Version1.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Version1.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /tmp/tmp.v25c7svKbI/main.cpp > CMakeFiles/Version1.dir/main.cpp.i

CMakeFiles/Version1.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Version1.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /tmp/tmp.v25c7svKbI/main.cpp -o CMakeFiles/Version1.dir/main.cpp.s

# Object files for target Version1
Version1_OBJECTS = \
"CMakeFiles/Version1.dir/main.cpp.o"

# External object files for target Version1
Version1_EXTERNAL_OBJECTS =

Version1: CMakeFiles/Version1.dir/main.cpp.o
Version1: CMakeFiles/Version1.dir/build.make
Version1: CMakeFiles/Version1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/tmp/tmp.v25c7svKbI/cmake-build-debug-remote-host-1/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Version1"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Version1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Version1.dir/build: Version1

.PHONY : CMakeFiles/Version1.dir/build

CMakeFiles/Version1.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Version1.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Version1.dir/clean

CMakeFiles/Version1.dir/depend:
	cd /tmp/tmp.v25c7svKbI/cmake-build-debug-remote-host-1 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /tmp/tmp.v25c7svKbI /tmp/tmp.v25c7svKbI /tmp/tmp.v25c7svKbI/cmake-build-debug-remote-host-1 /tmp/tmp.v25c7svKbI/cmake-build-debug-remote-host-1 /tmp/tmp.v25c7svKbI/cmake-build-debug-remote-host-1/CMakeFiles/Version1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Version1.dir/depend

