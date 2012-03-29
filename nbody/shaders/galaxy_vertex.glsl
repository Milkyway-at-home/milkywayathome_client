
#version 150

uniform mat4 cameraToClipMatrix;
uniform mat4 modelToCameraMatrix;

in vec3 position;
out vec2 galaxyPoint;

void main()
{
    gl_Position = cameraToClipMatrix * (modelToCameraMatrix * vec4(position, 1.0f));
    galaxyPoint = position.xy;
}

