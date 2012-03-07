
#version 150

out vec4 outputColor;

uniform sampler2D textTexture;

in vec2 texCoord;

void main()
{
    outputColor = texture(textTexture, texCoord);
}

