set(MORPH "${CMAKE_SOURCE_DIR}/data/morph.asc")

foreach(t distance_transform params regression compartments stiffness clustering)
  add_executable(test_${t} test_${t}.cpp)
  target_link_libraries(test_${t} PRIVATE occupancy::core)
  target_include_directories(test_${t} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()

# unit tests, no morphology needed
add_test(NAME distance_transform COMMAND test_distance_transform)
add_test(NAME params COMMAND test_params)
add_test(NAME stiffness COMMAND test_stiffness)
add_test(NAME clustering COMMAND test_clustering)

# these read the real morphologies
add_test(NAME regression COMMAND test_regression ${MORPH})
add_test(NAME compartments COMMAND test_compartments ${MORPH})
