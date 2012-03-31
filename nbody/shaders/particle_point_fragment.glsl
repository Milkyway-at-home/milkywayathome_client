
#version 150

flat in vec4 color;
out vec4 outputColor;

void main()
{
    outputColor = vec4(color.rgb, 1.0f);
}

