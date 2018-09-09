project(OpenFab)
cmake_minimum_required(VERSION 3.11)

include(FindProtobuf)
find_package(Protobuf REQUIRED)
include_directories(${PROTOBUF_INCLUDE_DIR})

##TODO: move to global CMake?
#list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../ExternalLibs/libigl/shared/cmake)

set(OpenFab_Include_Dir ${PROJECT_SOURCE_DIR}/include)
set(OpenFab_Source_Dir ${PROJECT_SOURCE_DIR}/src)

file(GLOB_RECURSE OpenFab_Header ${OpenFab_Include_Dir}/*.hpp)
file(GLOB_RECURSE OpenFab_Source ${OpenFab_Source_Dir}/*.cpp)

#Main setup
add_executable(OpenFab ${OpenFab_Header} ${OpenFab_Source})
target_include_directories(OpenFab PRIVATE ${OpenFab_Include_Dir})

#Add external libraries
include(../cmake_utils/cmake_utils.txt)
add_includes(OpenFab)

#Add other projects
target_include_directories(OpenFab PRIVATE ../DataStructure/include)
target_include_directories(OpenFab PRIVATE ../Geometry/include)
target_include_directories(OpenFab PRIVATE ../AutodiffEngine/include)
target_include_directories(OpenFab PRIVATE ../Mesh/include)
target_include_directories(OpenFab PRIVATE ../NetworkCommunication/include)
target_include_directories(OpenFab PRIVATE ../FabTranslation/include)


#target_include_directories(OpenFab PRIVATE ${PROJECT_SOURCE_DIR}/../ExternalLibs/eigen-git-mirror/Eigen)
#target_include_directories(OpenFab PRIVATE ${PROJECT_SOURCE_DIR}/../ExternalLibs/eigen-git-mirror/unsupported)
#target_include_directories(OpenFab PRIVATE ${PROJECT_SOURCE_DIR}/../ExternalLibs/CppAD)
#target_include_directories(OpenFab PRIVATE ${PROJECT_SOURCE_DIR}/../ExternalLibs/CppADCodeGen/include)


set_target_properties(OpenFab
  PROPERTIES
  DEBUG_POSTFIX "d_${LIB_EXTENSION}"
  RELEASE_POSTFIX "_${LIB_EXTENSION}"
  LINKER_LANGUAGE CXX
)
set(CMAKE_CXX_STANDARD_LIBRARIES -ldl)
set(CMAKE_VERBOSE_MAKEFILE ON)

#target_link_libraries(OpenFab igl::core)
target_link_libraries(OpenFab
        glfw
        glad
        NetworkCommunication
        ${PROTOBUF_LIBRARY})
