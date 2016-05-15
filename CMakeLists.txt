cmake_minimum_required(VERSION 3.0)
project(deflate)

if (MSVC)
    set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Supported build configurations" FORCE)

    # Reset flags

    set(MSVC_FLAGS_COMMON "/Zi /Zo /W4") # Always include debug information and set warning level 4

    set(CMAKE_C_FLAGS   "${MSVC_FLAGS_COMMON}")
    set(CMAKE_CXX_FLAGS "${MSVC_FLAGS_COMMON} /EHsc")

    set(MSVC_FLAGS_DEBUG   "/MTd /Od /D_DEBUG")
    set(MSVC_FLAGS_RELEASE "/MT /Ox /Ob2 /Oi /Ot /DNDEBUG /FAs")

    set(CMAKE_C_FLAGS_DEBUG     "${MSVC_FLAGS_DEBUG}")
    set(CMAKE_C_FLAGS_RELEASE   "${MSVC_FLAGS_RELEASE}")
    set(CMAKE_CXX_FLAGS_DEBUG   "${MSVC_FLAGS_DEBUG}")
    set(CMAKE_CXX_FLAGS_RELEASE "${MSVC_FLAGS_RELEASE}")

    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /DEBUG")
    set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} /ignore:4221") # warning LNK4221: This object file does not define any previously undefined public symbols, so it will not be used by any link operation that consumes this library

else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -pedantic -std=c++1z -g")
    if (MINGW)
        include_directories(${CMAKE_CURRENT_SOURCE_DIR}/third_party/mingw)
    endif()
    #TODO: Make sure boost filesystem is available
endif()

if (WIN32)
    add_definitions("-DWIN32 -D_SCL_SECURE_NO_WARNINGS -DUNICODE -D_UNICODE")
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    # On windows use: cmake .. -TLLVM-vs2014 -G"Visual Studio 14 2015 Win64"
    add_definitions("-Wno-missing-braces -Wall -Wextra")
endif()

enable_testing()
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --force-new-ctest-process --output-on-failure --build-config "$<CONFIGURATION>")

add_subdirectory(src)
include_directories(src)

add_executable(core_tests core_tests.cpp)
target_link_libraries(core_tests deflate_core)
add_test(core_tests core_tests)

add_executable(deflate main.cpp)
target_link_libraries(deflate deflate_core)