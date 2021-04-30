
#version 150

uniform mat4 cameraToClipMatrix;
uniform mat4 modelToCameraMatrix;

in vec4 position;
in vec4 inputColor;

flat out vec4 color;

void main()
{
    gl_Position = cameraToClipMatrix * modelToCameraMatrix * position;
    color = inputColor;
}

