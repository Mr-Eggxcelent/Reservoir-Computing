#pragma once
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "camera.h"
#include <iostream>


// //great solution by N0vember https://stackoverflow.com/questions/7676971/pointing-to-a-function-that-is-a-class-member-glfw-setkeycallback
// The data structure was adapted from Hazel Engine by Yan Chernikov https://github.com/TheCherno/Hazel/blob/master/Hazel/src/Platform/Windows/WindowsWindow.cpp
class Window
{

    GLFWwindow* _window;

    struct WindowData
    {
        std::string _Title;
        unsigned int _Width;
        unsigned int _Height;

        Camera& _camera;
        float _lastX = (float)_Width / 2.0;
        float _lastY = (float)_Height / 2.0;
        bool _firstMouse = true;

    public:


        WindowData(Camera& camera, unsigned int Width, unsigned int Height)
            :_camera(camera), _Width(Width),_Height(Height)
        {

        }
        // glfw: whenever the window size changed (by OS or user resize) this callback function executes
        // ---------------------------------------------------------------------------------------------
        void framebuffer_size_callback(GLFWwindow* window, int width, int height)
        {
            // make sure the viewport matches the new window dimensions; note that width and 
            // height will be significantly larger than specified on retina displays.
            glViewport(0, 0, width, height);
        }

        // glfw: whenever the mouse moves, this callback is called
        // -------------------------------------------------------
        void mouse_callback(GLFWwindow* window, double xpos, double ypos)
        {
            if (_firstMouse)
            {
                _lastX = xpos;
                _lastY = ypos;
                _firstMouse = false;
            }

            float xoffset = xpos - _lastX;
            float yoffset = _lastY - ypos; // reversed since y-coordinates go from bottom to top

            _lastX = xpos;
            _lastY = ypos;

            _camera.ProcessMouseMovement(xoffset, yoffset);
        }

        // glfw: whenever the mouse scroll wheel scrolls, this callback is called
        // ----------------------------------------------------------------------
        void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
        {
            _camera.ProcessMouseScroll(yoffset);
        }

    };

    WindowData m_Data;

public:

    Window(Camera& camera,unsigned int Width, unsigned int Height)
        :m_Data(camera,Width,Height)
    {

    }

    void initWindow()
    {

        // glfw: initialize and configure
       // ------------------------------
        glfwInit();
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

        _window = glfwCreateWindow(m_Data._Width, m_Data._Height, "LearnOpenGL", NULL, NULL);


        // glfw window creation
        // --------------------
        if (_window == NULL)
        {
            std::cout << "Failed to create GLFW window" << std::endl;
            glfwTerminate();
        }
        

        glfwSetWindowUserPointer(_window, &m_Data);

        auto frambuffer_func = [](GLFWwindow* w, int width, int height)
        { 
            static_cast<WindowData*>(glfwGetWindowUserPointer(w))->framebuffer_size_callback(w, width, height);
        };

        auto mouse_func = [](GLFWwindow* w, double xpos, double ypos)
        {
            static_cast<WindowData*>(glfwGetWindowUserPointer(w))->mouse_callback(w, xpos, ypos);
        };

        auto scroll_func = [](GLFWwindow* w, double xoffset, double yoffset)
        {
            static_cast<WindowData*>(glfwGetWindowUserPointer(w))->scroll_callback(w, xoffset, yoffset);
        };

        glfwMakeContextCurrent(_window);
        glfwSetFramebufferSizeCallback(_window, frambuffer_func);
        glfwSetCursorPosCallback(_window, mouse_func);
        glfwSetScrollCallback(_window, scroll_func);

        // tell GLFW to capture our mouse
        glfwSetInputMode(_window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);


    }

    GLFWwindow* getHandle()
    {
        return _window;
    }

    unsigned int getWidth() const { return m_Data._Width; }
    unsigned int getHeight() const { return m_Data._Height; }

};




