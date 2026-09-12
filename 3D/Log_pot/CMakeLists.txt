add_executable(LogPotExec source/LogPot_3D.cpp source/Initial_conditions.cpp)
    
target_include_directories(LogPotExec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../VecUtils/include ${CMAKE_CURRENT_SOURCE_DIR}/include)
    
target_link_libraries(LogPotExec PRIVATE VecUtilsLib)

add_executable(InitialConditionsExec source/Initial_conditions.cpp)
    
target_compile_definitions(InitialConditionsExec PRIVATE BUILD_INITIAL_CONDITIONS_MAIN)

target_include_directories(InitialConditionsExec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../VecUtils/include ${CMAKE_CURRENT_SOURCE_DIR}/include)
    
target_link_libraries(InitialConditionsExec PRIVATE VecUtilsLib)