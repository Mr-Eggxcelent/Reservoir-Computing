#version 330 core

layout (location = 0) in vec3 vertexPosition;
layout (location = 1) in vec3 vertexNormal;
layout (location = 2) in vec2 vertexTexCoord;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

// varyings (output)
out vec3 esVertex;
out vec3 esNormal;
out vec2 texCoord0;

void main()
{
    texCoord0 = vertexTexCoord;
    esNormal = mat3(transpose(inverse(model))) * vertexNormal;  
    esVertex = vec3(model * vec4(vertexPosition, 1.0));
    gl_Position =  projection * view * model * vec4(vertexPosition, 1.0);
}