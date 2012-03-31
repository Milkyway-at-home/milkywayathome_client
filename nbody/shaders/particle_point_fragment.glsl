
#version 150

flat in vec4 color;
out vec4 outputColor;

void main()
{
    // darken the ignored particles
    vec4 modColor = (color.a == 0.0f) ? (0.5f * color) : color;
    modColor.a = 1.0f;

    outputColor = modColor;
}

