
#version 150

out vec4 outputColor;

uniform sampler2D textTexture;

in vec2 texCoord;

void main()
{
    // swizzle red component to alpha since GL_ALPHA removed
    outputColor = texture(textTexture, texCoord).xxxx;
}

