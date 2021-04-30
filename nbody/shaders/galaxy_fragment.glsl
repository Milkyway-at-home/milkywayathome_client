
#version 150

uniform sampler2D galaxyTexture;
uniform float invGalaxyDiameter;

out vec4 outputColor;

in vec2 galaxyPoint;

/*
const vec4 white = vec4(1.0f, 1.0f, 1.0f, 1.0f);
const vec4 red = vec4(1.0f, 0.0f, 0.0f, 1.0f);
const vec4 green = vec4(0.0f, 1.0f, 0.0f, 1.0f);
*/

void main()
{
    vec2 texCoord = vec2(invGalaxyDiameter) * galaxyPoint + vec2(0.5f);
    vec4 tex = texture(galaxyTexture, texCoord);

    float intensity = (1.0f / 3.0f) * (tex.x + tex.y + tex.z);
    tex.w += intensity;

    outputColor = tex;

}

