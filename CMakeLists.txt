# ========== 基础配置 =========
cmake_minimum_required(VERSION 3.16)
project(PoissonSolvers 
        VERSION 1.2.0 
        LANGUAGES CXX
        DESCRIPTION "Finite Element Method for 2D Poisson Equation")

# 设置 C++17 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # PIC属性

# 编译器警告设置
if(MSVC)
    add_compile_options(/W4 /WX /permissive-)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
else()
    add_compile_options(-Wall -pedantic -fdiagnostics-color=always)
    add_compile_options(-march=native -O3)
endif()

# 设置输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ======= 核心库 (重量库) ======
add_library(PoissonCore SHARED)

# Windows 平台动态库导出控制
if(WIN32 AND BUILD_SHARED_LIBS)
    target_compile_definitions(PoissonCore PRIVATE POISSONCORE_EXPORTS)
endif()

# 使用 target_sources 添加源文件
target_sources(PoissonCore PRIVATE
    src/ReadInput.cpp
    src/ReadInput.h
    src/Timer.cpp
    src/Timer.h
    src/Help.cpp
    src/Help.h
    src/Mesh.cpp
    src/Mesh.h
    src/ExpressionEvaluator.cpp
    src/ExpressionEvaluator.h
    src/ElementCalculator.cpp
    src/ElementCalculator.h
    src/PoissonSolverBase.cpp
    src/PoissonSolverBase.h
    src/PoissonSolverDense.cpp
    src/PoissonSolverDense.h
    src/PoissonSolverSparse.cpp
    src/PoissonSolverSparse.h
    src/Core_Export.h
)

# 设置头文件目录
target_include_directories(PoissonCore PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/lib>
    $<INSTALL_INTERFACE:include>
)

# ============= 客户端可执行文件 (轻量客户端) ===============
add_executable(PoissonClient)
target_sources( PoissonClient PRIVATE src/main.cpp ) # 添加源文件
target_link_libraries(PoissonClient PRIVATE PoissonCore) # 链接核心库

# ============== 单元测试可执行文件 ============
enable_testing()

# 安装 GTest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
FetchContent_MakeAvailable(googletest)

# 添加测试
file(GLOB TEST_SOURCES CONFIGURE_DEPENDS 
    ./test/test_ExpressionEvaluator.cpp  
    ./test/test_Mesh.cpp  
    ./test/test_ReadInput.cpp
    )
add_executable(PoissonTests ${TEST_SOURCES})
target_link_libraries(PoissonTests PRIVATE PoissonCore gtest_main)
add_test(NAME AllTests COMMAND PoissonTests)

# ============ 依赖管理 ============
# 查找依赖库
find_package(Eigen3 3.3 REQUIRED)
find_package(Boost 1.71 REQUIRED)
find_package(GTest REQUIRED)

# 链接依赖
target_link_libraries(PoissonCore PUBLIC
    Eigen3::Eigen
    Boost::headers
)

# ========== 实用目标 ===============
# 客户端运行目标
add_custom_target(run
    COMMAND $<TARGET_FILE:PoissonClient>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Running Poisson Solver"
)

# 测试运行目标
add_custom_target(runtests
    COMMAND $<TARGET_FILE:PoissonTests>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Running Unit Tests"
)

# =========== 安装设置 ================
# 安装目标
install(TARGETS PoissonClient EXPORT PoissonSolversConfig
    RUNTIME DESTINATION bin
    COMPONENT applications
)

install(TARGETS PoissonCore EXPORT PoissonSolversConfig
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    COMPONENT libraries
)

install(DIRECTORY src/ 
    DESTINATION include/PoissonSolver
    COMPONENT headers
    FILES_MATCHING 
        PATTERN "*.h"
        PATTERN "*.hpp"
)

# 配置文件
install(EXPORT PoissonSolversConfig
    FILE PoissonSolversConfig.cmake
    NAMESPACE PoissonSolvers::
    DESTINATION lib/cmake/PoissonSolvers
    COMPONENT development
)
