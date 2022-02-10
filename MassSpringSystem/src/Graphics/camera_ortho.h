#pragma once

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <vector>

// Defines several possible options for camera movement. Used as abstraction to stay away from window-system specific input methods
enum Camera_Movement {
	ZOOM_IN,
	ZOOM_OUT,
	UP,
	DOWN,
	LEFT,
	RIGHT
};

// Default camera values
const float SPEED = 10.0f;
const float SENSITIVITY = 0.025f;
const float ZOOM = 1.0f;


// An abstract camera class that processes input and calculates the corresponding Euler Angles, Vectors and Matrices for use in OpenGL
class Camera
{

public:
	// Constructor with vectors
	Camera(float const left,
		   float const right,
		   float const bottom,
		   float const top,
		   glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f)) 
		: MovementSpeed(SPEED), MouseSensitivity(SENSITIVITY), Zoom(ZOOM)
		,_left(left),_right(right),_bottom(bottom),_top(top), Position(position)
	{
		_projection = glm::ortho(_left, _right / Zoom, _bottom, _top / Zoom, -1.0f, 1.0f);
	}

	// Constructor with scalar values
	Camera(float const left,
		   float const right,
		   float const bottom,
		   float const top,
		   float posX=0,float posY=0,float posZ=0) 
		: MovementSpeed(SPEED), MouseSensitivity(SENSITIVITY), Zoom(ZOOM),
		 _left(left), _right(right), _bottom(bottom), _top(top), Position(glm::vec3(posX, posY, posZ))
	{

		_projection = glm::ortho(_left, _right / Zoom, _bottom, _top / Zoom, -1.0f, 1.0f);
	}

	// Returns the view matrix calculated using Euler Angles and the LookAt Matrix
	glm::mat4 GetViewMatrix()
	{
		glm::mat4 transform = glm::translate(glm::mat4(1.0f), Position) *
			glm::rotate(glm::mat4(1.0f), rotation, glm::vec3(0, 0, 1));


		return glm::inverse(transform);
	}


	// Processes input received from any keyboard-like input system. Accepts input parameter in the form of camera defined ENUM (to abstract it from windowing systems)
	void ProcessKeyboard(Camera_Movement direction,float deltaTime)
	{
		float velocity = MovementSpeed * deltaTime;

		if (direction == UP)
			Position.y += velocity/ Zoom;
		if (direction == DOWN)
			Position.y -=  velocity/ Zoom;
		if (direction == LEFT)
			Position.x -= velocity/ Zoom;
		if (direction == RIGHT)
			Position.x += velocity/ Zoom;
	}


	//Refer to this https://stackoverflow.com/questions/67906914/opengl-move-2d-orthographic-camera-with-mouse
	//and this https://www.youtube.com/watch?v=ZQ8qtAizis4 by Javidx9
	void ProcessMouseMovement(double width, double height, double xoffset, double yoffset)
	{
		if (mouse_held) {

			xoffset *= MouseSensitivity;
			yoffset *= MouseSensitivity;

			Position.x -= ((float)xoffset) /(Zoom);
			Position.y -= ((float)yoffset) /(Zoom);

		}
	}

	// Processes input received from a mouse scroll-wheel event. Only requires input on the vertical wheel-axis
	//https://stackoverflow.com/questions/21561724/opengl-google-maps-style-2d-camera-zoom-to-mouse-cursor
	// This is a practical implementation https://stackoverflow.com/questions/21561724/opengl-google-maps-style-2d-camera-zoom-to-mouse-cursor by genpfault
    //This is the theory https://www.youtube.com/watch?v=ZQ8qtAizis4    Javidx9

	void ProcessMouseScroll(const double& xpos,const double& ypos,const double& yoffset)
	{

		glm::ivec4 viewport;
		glGetIntegerv(GL_VIEWPORT, &viewport[0]);

		glm::vec3 fWorldBefore = glm::unProject(
			glm::vec3(xpos, ypos, 0.0),
			GetViewMatrix(),
			_projection,
			viewport);

		if (Zoom >= 1.0f && Zoom <= 45.0f)
			Zoom += yoffset;
		if (Zoom <= 1.0f)
			Zoom = 1.0f;
		if (Zoom >= 45.0f)
			Zoom = 45.0f;

		_projection = glm::ortho(_left, _right/ Zoom, _bottom, _top / Zoom, -1.0f, 1.0f);

		glm::vec3 fWorldAfter = glm::unProject(
			glm::vec3(xpos, ypos, 0.0),
			GetViewMatrix(),
			_projection,
			viewport);

		Position.x += (fWorldBefore[0] - fWorldAfter[0]);
		Position.y += (fWorldBefore[1] - fWorldAfter[1]);
	}


	void setProjectionMatrix(glm::mat4& projection)
	{
		_projection = projection;
	}

	glm::mat4 getProjectionMatrix()
	{
		return _projection;
	}


private:
	glm::mat4 _projection;

	float const _left;
	float const _right;
	float const _bottom;
	float const _top;


public:
	// Camera Attributes
	glm::vec3 Position;

	bool mouse_held = false;

	// Camera options
	float MovementSpeed;
	float MouseSensitivity;
	float Zoom;
	float rotation = 0.0f;

};