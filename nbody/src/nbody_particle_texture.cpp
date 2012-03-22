/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

#include "nbody_gl_includes.h"
#include "nbody_particle_texture.h"
#include <cmath>


/**
* EvalHermite(float pA, float pB, float vA, float vB, float u)
* @brief Evaluates Hermite basis functions for the specified coefficients.
*/
static float evalHermite(float pA, float pB, float vA, float vB, float u)
{
    float u2 = u * u;
    float u3 = u2 * u;
    float B0 = 2.0f * u3 - 3.0f * u2 + 1.0f;
    float B1 = -2.0f * u3 + 3.0f * u2;
    float B2 = u3 - 2.0f * u2 + u;
    float B3 = u3 - u;
    return B0 * pA + B1 * pB + B2 * vA + B3 * vB;
}

/*
static float reducedEvalHermite(float u)
{
    return 2.0f * (u * u * u) - 3.0f * (u * u) + 1.0f;
}
*/

// from Nvidia OpenCL examples
static unsigned char* createGaussianMap(int n)
{
    float incr = 2.0f / (float) n;
    int i = 0;
    int j = 0;
    float y = -1.0f;
    float* m = new float[2 * n * n];
    unsigned char* b = new unsigned char[4 * n * n];

    for (int yi = 0; yi < n; ++yi, y += incr)
    {
        float y2 = y * y;
        float x = -1.0f;
        for (int xi = 0; xi < n; ++xi, x += incr, i += 2, j += 4)
        {
            float dist = sqrtf(x * x + y2);
            if (dist > 1.0f)
                dist = 1.0f;

            //m[i + 1] = m[i] = reducedEvalHermite(dist);
            //b[j] = (unsigned char)(255.0f * m[i]);
            m[i + 1] = m[i] = evalHermite(1.0f, 0.0f, 0.0f, 0.0f, dist);
            b[j] = b[j + 1] = b[j + 2] = b[j + 3] = (unsigned char) (255.0f * m[i]);
        }
    }

    delete[] m;

    return b;
}

GLuint createParticleTexture(int resolution)
{
    GLuint texture;
    unsigned char* data = createGaussianMap(resolution);

    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, resolution, resolution, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
    glGenerateMipmap(GL_TEXTURE_2D);

    delete[] data;

    return texture;
}

