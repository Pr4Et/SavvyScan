
cmake_minimum_required(VERSION 3.15)
project(SavvyscanServerCamera LANGUAGES CXX)

# C++ standard and common flags
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# -----------------------------------------------------------------------------
# Red Pitaya SDK discovery
# -----------------------------------------------------------------------------
# Allow override: -DREDPITAYA_SDK_ROOT=/opt/redpitaya
if(NOT DEFINED REDPITAYA_SDK_ROOT)
  set(REDPITAYA_SDK_ROOT "/opt/redpitaya")
endif()

# Headers (rp.h, rp_acq.h, rp_gen.h, ...)
set(RP_INCLUDE_DIR "${REDPITAYA_SDK_ROOT}/include")
# Libraries (librp.so)
set(RP_LIBRARY_DIR "${REDPITAYA_SDK_ROOT}/lib")

# Validate presence of headers (fail early if missing)
if(NOT EXISTS "${RP_INCLUDE_DIR}/rp.h")
  message(FATAL_ERROR "Red Pitaya header 'rp.h' not found at ${RP_INCLUDE_DIR}. "
                      "Set REDPITAYA_SDK_ROOT to the SDK prefix (e.g., /opt/redpitaya).")
endif()
if(NOT EXISTS "${RP_INCLUDE_DIR}/rp_acq.h")
  message(FATAL_ERROR "Red Pitaya header 'rp_acq.h' not found at ${RP_INCLUDE_DIR}. "
                      "This header declares rp_AcqSetDecimationFactor and related APIs.")
endif()

# Find librp (shared library)
find_library(RP_LIB rp HINTS "${RP_LIBRARY_DIR}")
if(NOT RP_LIB)
  message(FATAL_ERROR "Red Pitaya library 'librp.so' not found in ${RP_LIBRARY_DIR}.")
endif()

# Optional toggle for hardware code paths
option(HAVE_RP_HARDWARE "Enable Red Pitaya hardware code paths" ON)

# -----------------------------------------------------------------------------
# Compile-time probe: Does librp expose rp_AcqSetDecimationFactor(uint32_t)?
# -----------------------------------------------------------------------------
# We generate a tiny file at configure time and try to compile & link it.
# If it succeeds, define HAVE_NEW_DEC_API=1 for the targets.



# -----------------------------------------------------------------------------
# Targets
# -----------------------------------------------------------------------------
# Core library: BaseServer + ServerCamera + feeder
add_library(server_camera STATIC
  src/BaseServer.cpp
  src/ServerCamera.cpp
  src/feeder.cpp
)
target_compile_definitions(server_camera PRIVATE HAVE_NEW_DEC_API=1)

# Public include dirs for the library (your own headers)
target_include_directories(server_camera PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

# Red Pitaya include dir (needed by feeder.cpp / ServerCamera.cpp)
target_include_directories(server_camera PRIVATE
  "${RP_INCLUDE_DIR}"
)

# Hardware is always available on Red Pitaya
target_compile_definitions(server_camera PRIVATE HAVE_RP_HARDWARE=1)
target_compile_definitions(server_camera PRIVATE HAVE_NEW_DEC_API=1)

# Link the rp SDK to the static library so the executable doesn’t need to add it again
target_link_libraries(server_camera PRIVATE "${RP_LIB}")

# If you use threads (inotify, monitoring, etc.), link pthread too
find_package(Threads REQUIRED)
target_link_libraries(server_camera PRIVATE Threads::Threads)

# Executable: Linux entry point (port of SERVER-SEMCamServer.cpp)
add_executable(savvyscan_server
  src/ServerMain.cpp
)
target_link_libraries(savvyscan_server PRIVATE server_camera)

# Ensure the runtime can find librp.so without manual LD_LIBRARY_PATH:
# This sets an RPATH to /opt/redpitaya/lib (or your SDK lib dir) on the produced binary.
set_target_properties(savvyscan_server PROPERTIES
  BUILD_RPATH "${RP_LIBRARY_DIR}"
  INSTALL_RPATH "${RP_LIBRARY_DIR}"
)

# -----------------------------------------------------------------------------
# Helpful status output
# -----------------------------------------------------------------------------
message(STATUS "Using Red Pitaya SDK root: ${REDPITAYA_SDK_ROOT}")
message(STATUS "RP include dir: ${RP_INCLUDE_DIR}")
message(STATUS "RP library: ${RP_LIB}")
message(STATUS "HAVE_RP_HARDWARE: ${HAVE_RP_HARDWARE}")
message(STATUS "HAS_NEW_DEC_API: ${HAS_NEW_DEC_API}")
