# ===========================================================================
# EXERCISE 21 — CMake for a mixed Fortran + C + MPI + OpenMP project
# ===========================================================================
#
# GOAL
#   Build a small library that looks structurally like an ESM component:
#   Fortran modules with a dependency chain between them, a C kernel, MPI,
#   OpenMP, an installable target with exported CMake config, and tests.
#
# WHY (DKRZ)
#   "Scientific software design" is on the posting's wanted list, and in
#   practice a large part of an RSE's week is build systems. Specifically:
#
#     - Fortran MODULE dependencies are not textual includes. If a.f90 uses
#       module B, then b.f90 must be COMPILED first, because compiling it is
#       what produces b.mod. Make cannot infer that from #include scanning
#       the way it does for C. CMake's Fortran support does infer it — and
#       understanding what it is doing for you is what lets you fix it when
#       it goes wrong.
#     - Every model you are asked to port has a build system that predates
#       you, and "adapt complex models for European HPC systems" starts with
#       making the thing compile with a compiler it has never met.
#     - Mixed-language linking has a rule (link with the Fortran driver) and
#       an ABI (name mangling, which bind(C) exists to avoid).
#
# TASKS
#   TODO 1  project() with the right languages, and a minimum version
#   TODO 2  find_package for MPI and OpenMP; link them properly
#   TODO 3  the library target, with modules going somewhere sane
#   TODO 4  the executable, and getting the LINKER language right
#   TODO 5  install + export so another project can find_package() this
#   TODO 6  enable_testing + add_test, and answer the questions below
#
# ACCEPTANCE
#   - `make build` configures and compiles clean
#   - `make test` runs and passes
#   - `make deps` shows the Fortran module dependency order CMake worked out
#   - you can explain why physics.f90 must compile before driver.f90 without
#     either file naming the other in a build rule
#
# HINTS
#   - Modern CMake is target-based. Use target_link_libraries(tgt PUBLIC
#     MPI::MPI_Fortran OpenMP::OpenMP_Fortran) — never bare -lmpi and never
#     CMAKE_Fortran_FLAGS for a dependency's flags.
#   - Fortran_MODULE_DIRECTORY controls where .mod files land. Without it,
#     parallel builds of two targets in the same directory race over them,
#     and you get an intermittent build failure that looks like a compiler
#     bug and is not.
#   - LINKER_LANGUAGE Fortran on the executable, so the Fortran runtime
#     gets linked. A C driver linking Fortran objects fails on _gfortran_*
#     symbols, which is the same lesson as Exercise 04's Makefile.
#   - Test the OpenMP and MPI wiring, not just the compile: a test that runs
#     with 2 ranks and 2 threads catches a broken find_package immediately.
# ===========================================================================

# TODO 1: cmake_minimum_required and project().
# Pick a minimum you can defend -- 3.20 gives you decent Fortran support and
# is old enough to exist on cluster modules. Declare Fortran AND C.
cmake_minimum_required(VERSION 3.20)
project(minicomponent
        VERSION 0.1.0
        DESCRIPTION "A miniature ESM component, for build-system practice"
        LANGUAGES Fortran C)

# Default to a real build type. Without this, a bare `cmake ..` gives you
# NO optimisation flags at all, and people then report that "the model is
# slow" -- a genuinely common and embarrassing bug.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
  message(STATUS "No CMAKE_BUILD_TYPE set, defaulting to Release")
endif()

# ---------------------------------------------------------------------------
# TODO 2: find MPI and OpenMP.
#
#   find_package(MPI REQUIRED COMPONENTS Fortran C)
#   find_package(OpenMP REQUIRED COMPONENTS Fortran C)
#
# These define the imported targets MPI::MPI_Fortran, MPI::MPI_C,
# OpenMP::OpenMP_Fortran, OpenMP::OpenMP_C. Linking those propagates the
# include dirs, the compile flags AND the link flags automatically.
# ---------------------------------------------------------------------------
# TODO 2: uncomment and complete
# find_package(MPI REQUIRED COMPONENTS Fortran C)
# find_package(OpenMP REQUIRED COMPONENTS Fortran C)

# ---------------------------------------------------------------------------
# TODO 3: the library.
#
# Source files, and their real dependency order:
#     src/kinds.f90      module kinds          (no dependencies)
#     src/grid.f90       module grid_mod       uses kinds
#     src/physics.f90    module physics_mod    uses kinds, grid_mod
#     src/ckernel.c      C, called from physics via bind(C)
#
# List them in ANY order below -- CMake scans the `use` statements and works
# out the compile order itself. Deliberately list them backwards and confirm
# it still builds; that is the point of TODO 6a.
# ---------------------------------------------------------------------------
add_library(minicomponent
  src/physics.f90
  src/grid.f90
  src/kinds.f90
  src/ckernel.c
)

# TODO 3b: keep .mod files out of the source tree and out of other targets'
# way. Then tell consumers where to find them.
#
# set_target_properties(minicomponent PROPERTIES
#   Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
#
# target_include_directories(minicomponent PUBLIC
#   $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
#   $<INSTALL_INTERFACE:include>)

# TODO 3c: link the dependencies.
# target_link_libraries(minicomponent PUBLIC
#   MPI::MPI_Fortran MPI::MPI_C OpenMP::OpenMP_Fortran OpenMP::OpenMP_C)

# ---------------------------------------------------------------------------
# TODO 4: the executable.
#
# driver.f90 does `use mpi`, so it cannot compile until TODO 2 has wired up
# the MPI include path. The guard below keeps the shipped state building --
# the library compiles, proving the module-ordering lesson -- and the driver
# appears as soon as you complete TODO 2. Delete the guard once MPI is found.
# ---------------------------------------------------------------------------
if(TARGET MPI::MPI_Fortran)
  add_executable(driver src/driver.f90)
  target_link_libraries(driver PRIVATE minicomponent)

  # TODO 4b: force the Fortran linker driver. Try commenting this out once
  # the C file is doing real work and see what the link error looks like --
  # you want to recognise it instantly when it happens on a cluster at 6pm.
  # set_target_properties(driver PROPERTIES LINKER_LANGUAGE Fortran)
else()
  message(STATUS "")
  message(STATUS "  >> driver NOT built: MPI::MPI_Fortran does not exist.")
  message(STATUS "  >> Complete TODO 2 (find_package) to enable it.")
  message(STATUS "  >> The library still builds -- check `make deps` to see")
  message(STATUS "  >> the module order CMake derived on its own.")
endif()

# ---------------------------------------------------------------------------
# TODO 5: install and export.
#
# So that a DIFFERENT project can do find_package(minicomponent) and link
# against it. This is how ESM components get consumed by a coupled driver,
# and it is the difference between a library and a directory of files.
#
# include(GNUInstallDirs)
# install(TARGETS minicomponent driver
#         EXPORT minicomponentTargets
#         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
# install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/
#         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
# install(EXPORT minicomponentTargets
#         FILE minicomponentTargets.cmake
#         NAMESPACE minicomponent::
#         DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/minicomponent)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# TODO 6: tests. A build that is not tested is a build that is broken.
# ---------------------------------------------------------------------------
enable_testing()

# TODO 6: add tests that exercise the wiring, not just the compile.
#
# add_test(NAME runs_serial COMMAND driver)
#
# add_test(NAME runs_parallel
#          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
#                  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:driver> ${MPIEXEC_POSTFLAGS})
#
# set_tests_properties(runs_parallel PROPERTIES
#          ENVIRONMENT "OMP_NUM_THREADS=2"
#          PASS_REGULAR_EXPRESSION "ALL CHECKS PASSED")

message(STATUS "")
message(STATUS "  Exercise 21: work through the TODOs in this file.")
message(STATUS "  Verify with:  make build && make test && make deps")
message(STATUS "")

# ===========================================================================
# TODO 6 — write your answers in notes/day6.md
#
# (a) The source files are listed in reverse dependency order on purpose.
#     Explain how CMake still compiles them correctly. What does it do that
#     a hand-written Makefile does not? Now break it: make grid.f90 `use`
#     physics_mod as well, creating a cycle. What is the error, and why is
#     a circular module dependency unfixable rather than merely awkward?
#
# (b) Remove Fortran_MODULE_DIRECTORY and add a second library target in
#     this same directory that also compiles kinds.f90. Build with `make
#     -j8` repeatedly. What race appears, and how often? Why is an
#     intermittent build failure so much worse than a deterministic one?
#
# (c) Comment out LINKER_LANGUAGE Fortran. Paste the exact link error. Which
#     symbols are missing and where do they come from?
#
# (d) MPI::MPI_Fortran vs adding -I/opt/homebrew/include and -lmpi by hand.
#     List everything the imported target carries that the manual flags miss.
#     Why does that matter more on a cluster with modules than on a laptop?
#
# (e) A cluster gives you Intel compilers via `module load intel`. What
#     changes in this CMakeLists? (Answer: it should be nothing.) Now name
#     three things that commonly DO break in a real ESM build system when
#     the compiler changes, and how you would defend against each.
#
# (f) DKRZ context: ICON's build is autotools plus a large collection of
#     machine-specific config files, not CMake. Given the porting task in
#     the posting, what is the practical consequence for you as the new
#     person? Sketch what "adding support for a new HPC system" concretely
#     involves in such a build system.
# ===========================================================================
