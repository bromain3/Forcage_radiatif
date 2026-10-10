cmake_minimum_required(VERSION 4.1)  # Cmake minimal version for compatibility
project(code_C C)  # project name and languages

set(CMAKE_C_STANDARD 23)  # C version (2023)

# adding compile options for speed
if(MSVC)
    # Visual Studio
    add_compile_options(/O2 /arch:AVX2)
else()
    # GCC / Clang (Linux, macOS, MinGW on Windows)
    add_compile_options(-O3 -march=native -ffast-math -fopt-info-vec)  # we can deactivate fopt-info-vect if we don't want the optimisation info to appear during compilation
endif()


# Deactivating CURL stuff for improving compile time
set(BUILD_SHARED_LIBS        OFF CACHE BOOL "" FORCE)
set(BUILD_CURL_EXE           OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING            OFF CACHE BOOL "" FORCE)
set(CURL_DISABLE_TESTS       ON  CACHE BOOL "" FORCE)
set(CURL_DISABLE_LIBPSL     ON  CACHE BOOL "" FORCE)
set(CURL_DISABLE_LDAP       ON  CACHE BOOL "" FORCE)
set(CURL_DISABLE_LDAPS      ON  CACHE BOOL "" FORCE)
set(USE_NGHTTP2             OFF CACHE BOOL "" FORCE)
set(CURL_BROTLI             OFF CACHE BOOL "" FORCE)
set(CURL_ZSTD               OFF CACHE BOOL "" FORCE)
set(CURL_USE_LIBSSH2        OFF CACHE BOOL "" FORCE)
set(CURL_USE_LIBPSL         OFF CACHE BOOL "" FORCE)

# each OS has a builtin SSL system lib that we can use
if(WIN32)
    # Windows : SChannel is builtin
    set(CURL_USE_SCHANNEL   ON  CACHE BOOL "" FORCE)
    set(CURL_USE_OPENSSL    OFF CACHE BOOL "" FORCE)
elseif(APPLE)
    # macOS : Secure Transport is builtin
    set(CURL_USE_SECTRANSP  ON  CACHE BOOL "" FORCE)
    set(CURL_USE_OPENSSL    OFF CACHE BOOL "" FORCE)
else()
    # Linux : OpenSSL, available in each distribution (at least it should be)
    set(CURL_USE_OPENSSL    ON  CACHE BOOL "" FORCE)
endif()


# executable
add_subdirectory(curl)  # adding curl to the list of file to look in for libraries

add_executable(download dataDownloading/download.c dataDownloading/spe_ids_binding.h)  # defining a new executable and the files to compile for making it

target_link_libraries(download PRIVATE CURL::libcurl)  # linking curl libraries to the download code

add_executable(test_US_std_atm tests/test_US_std_atm.c athmosphericLayers/US_standard_atm/standardAtm.c)
target_compile_definitions(test_US_std_atm PRIVATE PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")  # allows to use the variable in the code so that the test code can call python easilly

add_executable(test_data_handler tests/test_data_handler.c dataHandler/dataHandler.c athmosphericLayers/US_standard_atm/standardAtm.c)
target_compile_definitions(test_data_handler PRIVATE PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")  # allows to use the variable in the code so that the test code can call python easilly

add_executable(test_integrator tests/test_integrator.c integrator/integrator.c athmosphericLayers/US_standard_atm/standardAtm.c dataHandler/dataHandler.c)
target_compile_definitions(test_integrator PRIVATE PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")  # allows to use the variable in the code so that the test code can call python easilly
