#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(
    ${CMAKE_SOURCE_DIR}
)

get_filename_component(DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${DIR_NAME})

add_executable(${DIR_NAME} "main.cpp")

#target_include_directories(${DIR_NAME} PRIVATE D:/soft/matlab/extern/include)
#target_link_directories(${DIR_NAME} PRIVATE D:/soft/matlab/extern/lib/win64/mingw64)
#target_link_libraries(${DIR_NAME} libmat libmx) #libmex libeng
