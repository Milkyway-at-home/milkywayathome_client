
#version 150

uniform sampler2D galaxyTexture;
uniform float invGalaxyDiameter;

out vec4 outputColor;

in vec2 galaxyPoint;

void main()
{
    vec2 scaledPoint = vec2(invGalaxyDiameter) * galaxyPoint;
    vec2 texCoord = scaledPoint + vec2(0.5f);
    vec4 tex = texture(galaxyTexture, texCoord);
    tex.w = 0.9f;
    outputColor = tex;
}

