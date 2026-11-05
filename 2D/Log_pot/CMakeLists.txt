add_executable(LogPotExec source/main.cpp source/LogPot.cpp source/SpiralPerturbation.cpp source/EffectivePotential.cpp source/Shu_DF.cpp)

target_include_directories(LogPotExec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../VecUtils/include ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(LogPotExec PRIVATE VecUtilsLib)



# Standalone Shu initial-condition generator (30000 stars, fixed seed 42).
add_executable(ShuDFExec source/Shu_DF.cpp)
target_compile_definitions(ShuDFExec PRIVATE BUILD_SHU_DF_MAIN)
target_include_directories(ShuDFExec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../VecUtils/include ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ShuDFExec PRIVATE VecUtilsLib)

# Optional visualization data export; reuses the existing equilibrium finder.
add_executable(LagrangeSnapshotsExec EXCLUDE_FROM_ALL
    source/LagrangeSnapshots.cpp source/EffectivePotential.cpp source/LogPot.cpp source/SpiralPerturbation.cpp)
target_include_directories(LagrangeSnapshotsExec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(LagrangeSnapshotsExec PRIVATE VecUtilsLib)
