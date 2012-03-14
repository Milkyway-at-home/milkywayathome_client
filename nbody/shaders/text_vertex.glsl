
#version 150

uniform mat4 cameraToClipMatrix;

in vec2 position;
in vec2 texPosition;

uniform sampler2D textTexture;
out vec2 texCoord;

void main()
{
    gl_Position = cameraToClipMatrix * vec4(position, 0.0f, 1.0f);
    texCoord = texPosition;
}

