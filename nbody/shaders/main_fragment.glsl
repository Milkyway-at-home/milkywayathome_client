
#version 150

//smooth in vec3 color;
flat in vec3 color;
out vec4 outputColor;

void main()
{
    //color = vec4(1.0f, 1.0f, 1.0f, 1.0f);
    //outputColor = vec4(Color, 1.0f);
    //outputColor = vec4(1.0f, 1.0f, 1.0f, 1.0f);
    /*
    float lerpValue = gl_FragCoord.y / 500.0f;

    outputColor = mix(vec4(1.0f, 1.0f, 1.0f, 1.0f),
                      vec4(0.2f, 0.2f, 0.2f, 1.0f),
                      lerpValue);
    */

    outputColor = vec4(color, 1.0f);
}

