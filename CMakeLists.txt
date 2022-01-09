cmake_minimum_required(VERSION "2.0.0")
project("MassSpringSystem")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHc /std:c++latest")
set(VENDOR_DIR "${CMAKE_CURRENT_SOURCE_DIR}/MassSpringSystem/vendor")
set(LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/MassSpringSystem/lib")

link_directories("${CMAKE_BINARY_DIR}/MassSpringSystem/lib")

math(EXPR PlatformArchitecture "8 * ${CMAKE_SIZEOF_VOID_P}")
set(SOURCES
	MassSpringSystem/src/Graphics/camera.h
	MassSpringSystem/src/Graphics/gnuplot-iostream.h
	MassSpringSystem/src/Graphics/render_text.h
	MassSpringSystem/src/Graphics/shader.h
	MassSpringSystem/src/Graphics/shapes.h
	MassSpringSystem/src/Graphics/sphere.cpp
	MassSpringSystem/src/Graphics/sphere.h
	MassSpringSystem/src/Graphics/window.h
	
	MassSpringSystem/src/Core/filereader.h
	MassSpringSystem/src/Core/Instrumentor.h
	MassSpringSystem/src/Core/main.cpp
	MassSpringSystem/src/Core/node.cpp
	MassSpringSystem/src/Core/node.h
	MassSpringSystem/src/Core/simulation.cpp
	MassSpringSystem/src/Core/simulation.h
	MassSpringSystem/src/Core/spring.cpp
	MassSpringSystem/src/Core/spring.h
	MassSpringSystem/src/Core/system.cpp
	MassSpringSystem/src/Core/system.h
	MassSpringSystem/src/Core/utility.cpp
	MassSpringSystem/src/Core/utility.h

)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/$<CONFIGURATION>-x${PlatformArchitecture}")
add_executable("${PROJECT_NAME}" ${SOURCES})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

#https://stackoverflow.com/questions/68014365/how-can-i-link-glfw-and-dear-imgui-in-cmake
#Erik Tellier

# OpenGL
find_package(OpenGL REQUIRED)

# GLM
find_package(GLM REQUIRED)
message(STATUS "GLM included at ${GLM_INCLUDE_DIR}")
target_include_directories(${PROJECT_NAME} PRIVATE "${GLM_INCLUDE_DIR}")

# GLFW
set(GLFW_DIR "${VENDOR_DIR}/GLFW")

set(GLFW_BUILD_EXAMPLES OFF CACHE INTERNAL "Build the GLFW example programs")
set(GLFW_BUILD_TESTS OFF CACHE INTERNAL "Build the GLFW test programs")
set(GLFW_BUILD_DOCS OFF CACHE INTERNAL "Build the GLFW documentation")
set(GLFW_INSTALL OFF CACHE INTERNAL "Generate installation target")

add_subdirectory("${GLFW_DIR}")
target_include_directories(${PROJECT_NAME} PRIVATE "${GLFW_DIR}/include")

#Freetype
set(FREETYPE_DIR "${VENDOR_DIR}/freetype")
add_subdirectory("${FREETYPE_DIR}")
target_include_directories(${PROJECT_NAME} PRIVATE "${FREETYPE_DIR}/include")

#Eigen
set(EIGEN_DIR "${VENDOR_DIR}/eigen")
target_include_directories(${PROJECT_NAME} PRIVATE "${EIGEN_DIR}")

#GLAD
set(GLAD_DIR "${VENDOR_DIR}/glad")
add_subdirectory("${GLAD_DIR}")
target_include_directories(${PROJECT_NAME} PUBLIC "${GLAD_DIR}/include")

#Delaunator
set(DELAUNATOR_DIR "${VENDOR_DIR}/delaunator-cpp")
target_include_directories(${PROJECT_NAME} PUBLIC "${DELAUNATOR_DIR}/include")

#Link Libraries
set(LIBS "-lg4tk" "glfw" "${GLFW_LIBRARIES}" "opengl32" "glad" "${CMAKE_DL_LIBS}" "${FREETYPE_LIBRARIES}" "freetype")
target_link_libraries(${PROJECT_NAME} ${LIBS})


# Define the link libraries
# MIT License
# Copyright (c) 2018 Tomasz Gałaj

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Create virtual folders to make it look nicer in VS
if(MSVC_IDE)
	# Macro to preserve source files hierarchy in the IDE
	macro(GroupSources curdir)
		file(GLOB children RELATIVE ${PROJECT_SOURCE_DIR}/${curdir} ${PROJECT_SOURCE_DIR}/${curdir}/*)

		foreach(child ${children})
			if(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/${curdir}/${child})
				GroupSources(${curdir}/${child})
			else()
				string(REPLACE "/" "\\" groupname ${curdir})
				string(REPLACE "src" "Sources" groupname ${groupname})
				source_group(${groupname} FILES ${PROJECT_SOURCE_DIR}/${curdir}/${child})
			endif()
		endforeach()
	endmacro()

	# Run macro
	GroupSources(src)
endif()


#https://stackoverflow.com/questions/66412965/how-to-set-the-output-and-intermediate-directory-in-visual-studio-c-projects-g?rq=1
#https://stackoverflow.com/questions/10671916/how-to-copy-dll-files-into-the-same-folder-as-the-executable-using-cmake
#https://shot511.github.io/2018-05-29-how-to-setup-opengl-project-with-cmake/
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
	add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_directory
		"${CMAKE_BINARY_DIR}/bin/$<CONFIGURATION>-x${PlatformArchitecture}"
		$<TARGET_FILE_DIR:${PROJECT_NAME}>)
endif()
		
		
		
		
		
		