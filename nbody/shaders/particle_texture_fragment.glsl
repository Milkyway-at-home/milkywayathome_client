
#version 150

uniform sampler2D particleTexture;

flat in vec4 color;
out vec4 outputColor;

void main()
{
    vec4 tex = texture(particleTexture, gl_PointCoord);
    outputColor = vec4(color.rgb, tex.a);
}

