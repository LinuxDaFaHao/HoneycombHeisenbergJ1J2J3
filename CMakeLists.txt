cmake_minimum_required(VERSION 3.14)
project(HoneycombHeisenberg)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
set(CMAKE_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src/")

# Build type.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

#if(NOT PHYSICAL_MODEL)
#    set(PHYSICAL_MODEL Fermion)
#endif()

if(${CMAKE_BUILD_TYPE} MATCHES Debug)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

option(QLTEN_TIMING_MODE "Timing mode." OFF)
if(QLTEN_TIMING_MODE)
    add_definitions(-DQLTEN_TIMING_MODE)
endif()

option(QLMPS_TIMING_MODE "Timing mode." OFF)
if(QLMPS_TIMING_MODE)
    add_definitions(-DQLMPS_TIMING_MODE)
endif()

option(QLTEN_MPI_TIMING_MODE "Timing mode." OFF)
if(QLTEN_MPI_TIMING_MODE)
    add_definitions(-DQLTEN_MPI_TIMING_MODE)
endif()

option(QLMPS_MPI_TIMING_MODE "Timing mode." OFF)
if(QLMPS_MPI_TIMING_MODE)
    add_definitions(-DQLMPS_MPI_TIMING_MODE)
endif()

option(COMPILE_FIX_CODE "compile the fixing code" OFF)

# Compilation and linking control.
set(CMAKE_CXX_STANDARD 17)
add_definitions(-Wall -g)

#add_compile_options(-O3)



if (APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
    if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(MATH_LIB_COMPILE_FLAGS "-I$ENV{MKLROOT}/include")
        #Need test
        set(MATH_LIB_LINK_FLAGS $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a $ENV{MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm -ldl)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
        # May not work
        set(MATH_LIB_LINK_FLAGS $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a $ENV{MKLROOT}/lib/libmkl_core.a -L$ENV{MKLROOT}/lib -L$ENV{CMPLR_ROOT}/mac/compiler/lib/ -liomp5 -lpthread -lm -ldl)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
        # Note as of Intel oneAPI 2021.2, "source /opt/intel/oneapi/compiler/latest/env/vars.sh"
        set(MATH_LIB_LINK_FLAGS $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a $ENV{MKLROOT}/lib/libmkl_core.a -L$ENV{MKLROOT}/lib -L$ENV{CMPLR_ROOT}/mac/compiler/lib/ -Wl, -rpath $ENV{CMPLR_ROOT}/mac/compiler/lib/libiomp5.dylib -liomp5 -lpthread -lm -ldl)
    endif()
elseif (APPLE) #ARM
    add_definitions(-DUSE_OPENBLAS)
    set(BLA_VENDOR OpenBLAS)
    #        FIND_PACKAGE(BLAS REQUIRED)
    #        FIND_PACKAGE(LAPACK REQUIRED)
    set(OpenBLAS_ROOT "/opt/homebrew/opt/openblas/")
    set(Lapack_ROOT "/opt/homebrew/opt/lapack")
    message(${OpenBLAS_ROOT})
    set(OpenBLAS_INCLUDE_DIRS "${OpenBLAS_ROOT}/include")
    set(OpenBLAS_LIBRARIES "${OpenBLAS_ROOT}/lib/libblas.dylib")
    message(${OpenBLAS_LIBRARIES})
    set(MATH_LIB_COMPILE_FLAGS -I${OpenBLAS_INCLUDE_DIRS} -pthread)
    set(MATH_LIB_LINK_FLAGS ${OpenBLAS_LIBRARIES} ${OpenBLAS_ROOT}/lib/liblapack.dylib -lm -lpthread -ldl -fopenmp -lclapack)

elseif(UNIX)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "~/.local/")
    message("CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")

    #UNIX include APPLE, but we except it here
    if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(MATH_LIB_COMPILE_FLAGS "-I$ENV{MKLROOT}/include")
        set(MATH_LIB_LINK_FLAGS -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
        # Link the Intel's OpenMP library to avoid performance issue when the library calls the MKL's gesdd function.
        set(MATH_LIB_LINK_FLAGS -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -L$ENV{MKLROOT}/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -lpthread -lirc -lsvml -Wl,--no-as-needed -ldl)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
        set(MATH_LIB_LINK_FLAGS -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -L$ENV{MKLROOT}/lib/intel64 -liomp5 -lpthread -lm -ldl)
    endif()
endif()



find_path(TENSOR_HEADER_PATH "QLTEN")
find_path(MPS2_HEADER_PATH "QLMPS")
find_package(hptt)
find_package(Boost
        1.76.0
        REQUIRED
        COMPONENTS serialization mpi)
find_package(MPI REQUIRED)
#set omp flag
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(OMP_FLAGS -qopenmp)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(OMP_FLAGS -fopenmp)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(OMP_FLAGS -fopenmp)
endif()

#find openmp?
add_compile_options( ${MATH_LIB_COMPILE_FLAGS} ${OMP_FLAGS})

include_directories(  ${MPS2_HEADER_PATH}
        ${TENSOR_HEADER_PATH}
        ${hptt_INCLUDE_DIR}
        ${Boost_INCLUDE_DIR}
        ${MPI_CXX_HEADER_DIR}
        )

link_libraries(
        ${hptt_LIBRARY}
        ${MATH_LIB_LINK_FLAGS}
        Boost::serialization
        Boost::mpi
        ${MPI_CXX_LINK_FLAGS}
        ${MPI_mpi_LIBRARY}
)

add_executable(mpogen src/mpogen.cpp src/operators.cpp)
add_executable(dmrg src/dmrg.cpp src/myutil.cpp)
add_executable(measure1 src/measure1.cpp src/myutil.cpp src/operators.cpp)
add_executable(measuresf src/measure_spin_structure_factor.cpp src/myutil.cpp src/operators.cpp)
add_executable(contract_demo src/contract_demo.cpp)

add_executable(fixenv src/fix_env.cpp src/myutil.cpp)
add_executable(fixmps src/fix_mps.cpp src/myutil.cpp)
add_executable(continue src/continue.cpp src/myutil.cpp)

add_executable(xtrg src/xtrg.cpp src/myutil.cpp)

add_executable(test1 XTRG_common/xy_model_test.cpp)
add_executable(test2 XTRG_common/transverse_ising_test.cpp)

add_executable(ising MonteCarlo_src/ising.cpp)
add_executable(ising_exchange MonteCarlo_src/ising_exchange.cpp)
add_executable(xy MonteCarlo_src/xymodel.cpp)
add_executable(heisenberg MonteCarlo_src/classical_heisenberg.cpp)
add_executable(statistic MonteCarlo_src/statistic.cpp)
add_executable(ising_statistic MonteCarlo_src/ising_statistic.cpp)
#add_executable(testani1 MonteCarlo_src/anisotropy_test1.cpp)
add_executable(test MonteCarlo_src/test.cpp)
add_executable(heisenberg_exchange MonteCarlo_src/heisenberg_exchange.cpp)
add_executable(statistic_hei_ex MonteCarlo_src/statistic_heisenberg_exMC.cpp)
add_executable(clock ClockModel/clock.cpp)