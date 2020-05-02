aux_source_directory(${CMAKE_SOURCE_DIR}/src FILE_SRC)

set(FILE_TEST
    ${FILE_SRC}
    test_main.cc
    test_finite.cc
    test_infinite.cc)

add_executable(unit_test
    ${FILE_TEST})

target_link_libraries(unit_test
    PRIVATE
    ${Boost_LIBRARIES}
    fmt-header-only)

target_include_directories(unit_test
    PRIVATE
    ${CMAKE_SOURCE_DIR}/include)