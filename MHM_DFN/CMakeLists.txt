file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(Mymain main.cpp ${sources} ${headers})

add_executable(MHM_DFN ${Mymain})

target_link_libraries(MHM_DFN pz)

