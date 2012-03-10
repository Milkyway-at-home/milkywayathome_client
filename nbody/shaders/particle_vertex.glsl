
#version 150

uniform mat4 cameraToClipMatrix;
uniform mat4 modelToCameraMatrix;

in vec3 position;
in vec3 inputColor;

uniform float pointSize;

flat out vec4 color;

void main()
{
    vec4 cameraPos = modelToCameraMatrix * vec4(position, 1.0f);
    gl_Position = cameraToClipMatrix * cameraPos;
    gl_PointSize = max(1.0f, pointSize / (1.0f - cameraPos.z));
    color = vec4(inputColor, 1.0f);
}

