cmake_minimum_required(VERSION 2.8)
project (SPOKEDARTS)

file(GLOB spokedarts_source
	"*.hpp"
	"*.cpp"
)

add_executable(spokedarts ${spokedarts_source})

add_definitions(-std=c++11) 

install(TARGETS spokedarts DESTINATION bin)
