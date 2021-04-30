
#version 150

uniform mat4 cameraToClipMatrix;

in vec4 position; // xy = space position, zw = st texture

uniform sampler2D textTexture;
out vec2 texCoord;

void main()
{
    gl_Position = cameraToClipMatrix * vec4(position.xy, 0.0f, 1.0f);
    texCoord = position.zw;
}

