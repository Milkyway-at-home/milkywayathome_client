
#version 150

uniform sampler2D particleTexture;

flat in vec4 color;
out vec4 outputColor;

void main()
{
    vec4 tex = texture(particleTexture, gl_PointCoord);

    // darken the ignored particles
    // TODO: Maybe a different color would be better? green?
    vec4 modColor = (color.a == 0.0f) ? (0.5f * color) : color;

    outputColor = vec4(modColor.rgb, tex.a);
}

