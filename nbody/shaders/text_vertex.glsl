
#version 150

uniform mat4 cameraToClipMatrix;

in vec3 position;
in vec2 texPosition;

uniform sampler2D textTexture;
out vec2 texCoord;

void main()
{
    gl_Position = cameraToClipMatrix * vec4(position, 1.0f);
    texCoord = texPosition;
}

