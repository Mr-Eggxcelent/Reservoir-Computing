#version 330 core

// Input vertex data, different for all executions of this shader.
layout(location = 0) in vec4 vertices;

// Output data ; will be interpolated for each fragment.
out vec2 TexCoords;

// Values that stay constant for the whole mesh.
uniform vec3 CameraRight_worldspace;
uniform vec3 CameraUp_worldspace;
uniform mat4 view;
uniform mat4 projection;
uniform vec3 center_pos;


void main()
{
	vec3 particleCenter_wordspace = center_pos;
	
	//vec3 vertexPosition_worldspace = center_pos;


	//vec3 vertexPosition_worldspace = 
		//particleCenter_wordspace
		//+ CameraRight_worldspace * vertices.x
		//+ CameraUp_worldspace * vertices.y;


	// Output position of the vertex
	gl_Position = projection * view * vec4(vertices.xy,0.0, 1.0);


	// UV of the vertex. No special space for this one.
	//TexCoords = vertices.xy + vec2(0.5, 0.5);
   TexCoords = vertices.zw;

}
