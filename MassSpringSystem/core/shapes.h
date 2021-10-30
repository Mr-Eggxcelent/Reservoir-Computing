#pragma once
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "shader.h"
#include "camera.h"

#include<array>
#include <iostream>
#include <memory>
#include "sphere.h"

//Can use polymorphism https://www.youtube.com/watch?v=kxKKHKSMGIg Javidx9
class SphereRender
{
    GLuint _vboId1;// IDs of VBO for vertex arrays
    GLuint _iboId1;// IDs of VBO for index array
    GLuint _vaoId1;
    
    std::unique_ptr<Shader>_shader;

    // sphere: min sector = 3, min stack = 2
    Sphere _sphere1;           // radius, sectors, stacks, smooth(default)
    glm::vec4 _material_ambient;

public:
    SphereRender()
        :_sphere1(0.08, 36, 18), _material_ambient(0.5f,0.5f,0.5f,1.0f)
    {
    }


    void initBuffer()
    {
        _shader = std::make_unique<Shader>("shader/sphere_shader.vs", "shader/sphere_shader.fs");
        // set attrib arrays using glVertexAttribPointer()
        int stride = _sphere1.getInterleavedStride();

        glGenBuffers(1, &_vboId1);
        glGenVertexArrays(1, &_vaoId1);
        glGenBuffers(1, &_iboId1);

        glBindVertexArray(_vaoId1);

        glBindBuffer(GL_ARRAY_BUFFER, _vboId1);
        glBufferData(GL_ARRAY_BUFFER, _sphere1.getInterleavedVertexSize(), _sphere1.getInterleavedVertices(), GL_STATIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _iboId1);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, _sphere1.getIndexSize(), _sphere1.getIndices(), GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, 0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, stride, (void*)(6 * sizeof(float)));
        glEnableVertexAttribArray(2);

    }


    void draw(glm::mat4& projection, glm::mat4& view, glm::vec3 position, glm::vec4 material_ambient)
    {
        _shader->use();
        glm::mat4 model = glm::mat4(1.0f);
        _shader->setMat4("projection", projection);
        _shader->setMat4("view", view);
        model = glm::translate(model, position);
        _shader->setMat4("model", model);

        _shader->setVec4("lightPosition", 1, 1, 1, 1);
        _shader->setVec4("lightAmbient", 0.5f, 0.5f, 0.5f, 1);
        _shader->setVec4("lightDiffuse", 0.7f, 0.7f, 0.7f, 1);
        _shader->setVec4("lightSpecular", 1.0f, 1.0f, 1.0f, 1);
        _material_ambient = material_ambient;
        _shader->setVec4("materialAmbient", _material_ambient);
        _shader->setVec4("materialDiffuse", 0.7f, 0.7f, 0.7f, 1);
        _shader->setVec4("materialSpecular", 0.4f, 0.4f, 0.4f, 1);
        _shader->setFloat("materialShininess", 30.0f);


        glBindVertexArray(_vaoId1);
        // draw center sphere
        glDrawElements(GL_TRIANGLES,            // primitive type
            _sphere1.getIndexCount(), // # of indices
            GL_UNSIGNED_INT,         // data type
            (void*)0);               // ptr to indices

    }

    void clearResources()
    {
        glDeleteBuffers(1, &_vboId1);
        glDeleteBuffers(1, &_iboId1);
        glDeleteBuffers(1, &_vaoId1);
        _vboId1 = _iboId1 = 0;

    }



};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//https://stackoverflow.com/questions/14486291/how-to-draw-line-in-opengl
//Modified Code originally written by stackoverflow user jackw11111
//Jan 7 '19 at 6:26
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
class LineRender
{
    std::unique_ptr<Shader>_shader;
    //std::vector<float> _vertices;
    std::array<float, 6>_vertices;
    GLuint _vbo;
    GLuint _vao;

    glm::vec3 lineColor;
    glm::vec3 _startPoint;
    glm::vec3 _endPoint;

    struct Vertex
    {
        float Position[3];
    };


public:

    LineRender() = default;

    //If you need to undestand this code watch this video
    //https://www.youtube.com/watch?v=5df3NvQNzUs by The Cherno
    void initBuffer()
    {
        _shader = std::make_unique<Shader>("shader/line_shader.vs", "shader/line_shader.fs");


        glGenVertexArrays(1, &_vao);
        glGenBuffers(1, &_vbo);

        glBindVertexArray(_vao);
        glBindBuffer(GL_ARRAY_BUFFER, _vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * 6, nullptr, GL_DYNAMIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const void*)offsetof(Vertex, Position));
        glEnableVertexAttribArray(0);
     

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
        

    }


    void draw(glm::mat4& projection, glm::mat4& view, glm::vec3 s_Point, glm::vec3 e_Point)
    {

        _startPoint = s_Point;
        _endPoint = e_Point;

        _vertices = {
          _startPoint.x, _startPoint.y, _startPoint.z,
          _endPoint.x, _endPoint.y, _endPoint.z
        };

       
        glBindBuffer(GL_ARRAY_BUFFER, _vbo);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(_vertices), _vertices.data());
       
        glm::mat4 model = glm::mat4(1.0f);


        _shader->use();
        _shader->setMat4("projection", projection);
        _shader->setMat4("model", model);
        _shader->setMat4("view", view);
        _shader->setVec3("color", lineColor);
        

        glBindVertexArray(_vao);
        glDrawArrays(GL_LINES, 0, 2);

    }


    void setColor(glm::vec3 color) {
        lineColor = color;
    }


    void clearResources()
    {
        glDeleteVertexArrays(1, &_vao);
        glDeleteBuffers(1, &_vbo);
    }


};
