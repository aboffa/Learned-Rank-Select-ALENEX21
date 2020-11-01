project(tests)

find_package(GTest REQUIRED)

add_executable(my_tests test.cpp)
target_link_libraries(my_tests ${GTEST_LIBRARIES} ${Boost_LIBRARIES} sdsl)