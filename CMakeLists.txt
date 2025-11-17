cmake_minimum_required(VERSION 3.10)
project(CPLEXExample)

set(CMAKE_CXX_STANDARD 11)

# 设置 CPLEX 路径
set(CPLEX_ROOT "/Applications/CPLEX_Studio1210")
set(CPLEX_INCLUDE_DIRS ${CPLEX_ROOT}/cplex/include ${CPLEX_ROOT}/concert/include)
set(CPLEX_LIBRARY_DIRS ${CPLEX_ROOT}/cplex/lib/x86-64_osx/static_pic ${CPLEX_ROOT}/concert/lib/x86-64_osx/static_pic)

# 添加自定义编译选项
add_compile_options(-DIL_STD)

# 添加可执行文件
add_executable(main main.cpp)
#
# 包含头文件目录
target_include_directories(main PUBLIC ${CPLEX_INCLUDE_DIRS})

# 链接库文件
target_link_libraries(main 
    ${CPLEX_LIBRARY_DIRS}/libilocplex.a
    ${CPLEX_LIBRARY_DIRS}/libcplex.a
    ${CPLEX_LIBRARY_DIRS}/libconcert.a
    pthread
    m
)

# 设置链接器搜索路径
set_target_properties(main PROPERTIES 
    LINK_FLAGS "-Wl,-rpath,${CPLEX_LIBRARY_DIRS}"
)