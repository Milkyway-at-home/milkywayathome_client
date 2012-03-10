
#version 150

uniform sampler2D particleTexture;

flat in vec4 color;
out vec4 outputColor;

in float pointSizeOut;

void main()
{
    //vec4 tex = texture(particleTexture, vec2(0.5f, 0.5f));
    //vec4 tex = texture(particleTexture, texPos);
    //outputColor = vec4(vec3(1.0f) - tex.xyz, tex.w);
    //outputColor = tex;
    //outputColor = color;

    vec4 tex = texture(particleTexture, gl_PointCoord);
    vec4 texs = (0.6f + 0.4f * color) * tex;

    //outputColor = texs * mix(vec4(0.0f, 0.2f, 0.2f, texs.w), vec4(0.2f, 0.7f, 0.7f, texs.w), texs.w);

    //outputColor = tex;
    outputColor = texs;
}

