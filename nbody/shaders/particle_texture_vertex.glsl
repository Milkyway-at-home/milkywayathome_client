
#version 150

//uniform sampler2D particleTexture;

uniform mat4 cameraToClipMatrix;
uniform mat4 modelToCameraMatrix;

in vec3 position;
in vec3 inputColor;

flat out vec4 color;

out float pointSizeOut;

void main()
{
    vec4 cameraPos = modelToCameraMatrix * vec4(position, 1.0f);
    gl_Position = cameraToClipMatrix * cameraPos;

    // gl_Point.size is gone?
    float pointSize = 100.0f;
    gl_PointSize = pointSizeOut = max(1.0f, pointSize / (1.0f - cameraPos.z));
    color = vec4(inputColor, 1.0f);
}

