add_executable(LogPotExec source/LogPot.cpp source/Dehnen_DF.cpp)

target_include_directories(LogPotExec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../VecUtils/include ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(LogPotExec PRIVATE VecUtilsLib)



# Standalone Shu initial-condition generator (30000 stars, fixed seed 42).
add_executable(ShuDFExec source/Dehnen_DF.cpp)
target_compile_definitions(ShuDFExec PRIVATE BUILD_SHU_DF_MAIN)
target_include_directories(ShuDFExec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../VecUtils/include ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ShuDFExec PRIVATE VecUtilsLib)
