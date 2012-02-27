
#version 150

//uniform mat4 viewMatrix;
//uniform mat4 projMatrix;

uniform mat4 cameraToClipMatrix;
uniform mat4 modelToCameraMatrix;

in vec3 position;
in vec3 inputColor;

//smooth out vec3 color;
flat out vec3 color;

void main()
{
    //gl_Position = projMatrix * viewMatrix * position;
    //gl_PointSize = 10.0f;
    //gl_Position = vec4(position, 1.0f);

    gl_Position = cameraToClipMatrix * (modelToCameraMatrix * vec4(position, 1.0f));
    color = inputColor;
}

