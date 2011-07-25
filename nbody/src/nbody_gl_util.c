/*
 * freeglut_font.c
 *
 * Bitmap and stroke fonts displaying.
 *
 * Copyright (c) 1999-2000 Pawel W. Olszta. All Rights Reserved.
 * Written by Pawel W. Olszta, <olszta@sourceforge.net>
 * Creation date: Thu Dec 16 1999
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * PAWEL W. OLSZTA BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "nbody_gl_util.h"
#include <stdlib.h>

typedef struct
{
    const char*     Name;         /* The source font name             */
    int             Quantity;     /* Number of chars in font          */
    int             Height;       /* Height of the characters         */
    const GLubyte** Characters;   /* The characters mapping           */

    float           xorig, yorig; /* Relative origin of the character */
} NBody_SFG_Font;

static const GLubyte Helvetica12_Character_000[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_001[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_002[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_003[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_004[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_005[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_006[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_007[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_008[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_009[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_010[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_011[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_012[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_013[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_014[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_015[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_016[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_017[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_018[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_019[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_020[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_021[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_022[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_023[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_024[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_025[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_026[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_027[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_028[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_029[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_030[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_031[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_032[] = {  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_033[] = {  3,  0,  0,  0,  0, 64,  0, 64, 64, 64, 64, 64, 64, 64,  0,  0,  0};
static const GLubyte Helvetica12_Character_034[] = {  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 80, 80, 80,  0,  0,  0};
static const GLubyte Helvetica12_Character_035[] = {  7,  0,  0,  0,  0, 80, 80, 80,252, 40,252, 40, 40,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_036[] = {  7,  0,  0,  0, 16, 56, 84, 84, 20, 56, 80, 84, 56, 16,  0,  0,  0};
static const GLubyte Helvetica12_Character_037[] = { 11,  0,  0,  0,  0,  0,  0,  0,  0, 17,128, 10, 64, 10, 64,  9,128,  4,  0, 52,  0, 74,  0, 74,  0, 49,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_038[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 57,  0, 70,  0, 66,  0, 69,  0, 40,  0, 24,  0, 36,  0, 36,  0, 24,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_039[] = {  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 64, 32, 96,  0,  0,  0};
static const GLubyte Helvetica12_Character_040[] = {  4,  0, 16, 32, 32, 64, 64, 64, 64, 64, 64, 32, 32, 16,  0,  0,  0};
static const GLubyte Helvetica12_Character_041[] = {  4,  0,128, 64, 64, 32, 32, 32, 32, 32, 32, 64, 64,128,  0,  0,  0};
static const GLubyte Helvetica12_Character_042[] = {  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 80, 32, 80,  0,  0,  0};
static const GLubyte Helvetica12_Character_043[] = {  7,  0,  0,  0,  0,  0, 16, 16,124, 16, 16,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_044[] = {  4,  0,  0, 64, 32, 32,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_045[] = {  8,  0,  0,  0,  0,  0,  0,  0,124,  0,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_046[] = {  3,  0,  0,  0,  0, 64,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_047[] = {  4,  0,  0,  0,  0,128,128, 64, 64, 64, 32, 32, 16, 16,  0,  0,  0};
static const GLubyte Helvetica12_Character_048[] = {  7,  0,  0,  0,  0, 56, 68, 68, 68, 68, 68, 68, 68, 56,  0,  0,  0};
static const GLubyte Helvetica12_Character_049[] = {  7,  0,  0,  0,  0, 16, 16, 16, 16, 16, 16, 16,112, 16,  0,  0,  0};
static const GLubyte Helvetica12_Character_050[] = {  7,  0,  0,  0,  0,124, 64, 64, 32, 16,  8,  4, 68, 56,  0,  0,  0};
static const GLubyte Helvetica12_Character_051[] = {  7,  0,  0,  0,  0, 56, 68, 68,  4,  4, 24,  4, 68, 56,  0,  0,  0};
static const GLubyte Helvetica12_Character_052[] = {  7,  0,  0,  0,  0,  8,  8,252,136, 72, 40, 40, 24,  8,  0,  0,  0};
static const GLubyte Helvetica12_Character_053[] = {  7,  0,  0,  0,  0, 56, 68, 68,  4,  4,120, 64, 64,124,  0,  0,  0};
static const GLubyte Helvetica12_Character_054[] = {  7,  0,  0,  0,  0, 56, 68, 68, 68,100, 88, 64, 68, 56,  0,  0,  0};
static const GLubyte Helvetica12_Character_055[] = {  7,  0,  0,  0,  0, 32, 32, 16, 16, 16,  8,  8,  4,124,  0,  0,  0};
static const GLubyte Helvetica12_Character_056[] = {  7,  0,  0,  0,  0, 56, 68, 68, 68, 68, 56, 68, 68, 56,  0,  0,  0};
static const GLubyte Helvetica12_Character_057[] = {  7,  0,  0,  0,  0, 56, 68,  4,  4, 60, 68, 68, 68, 56,  0,  0,  0};
static const GLubyte Helvetica12_Character_058[] = {  3,  0,  0,  0,  0, 64,  0,  0,  0,  0, 64,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_059[] = {  3,  0,  0,128, 64, 64,  0,  0,  0,  0, 64,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_060[] = {  7,  0,  0,  0,  0,  0, 12, 48,192, 48, 12,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_061[] = {  7,  0,  0,  0,  0,  0,  0,124,  0,124,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_062[] = {  7,  0,  0,  0,  0,  0, 96, 24,  6, 24, 96,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_063[] = {  7,  0,  0,  0,  0, 16,  0, 16, 16,  8,  8, 68, 68, 56,  0,  0,  0};
static const GLubyte Helvetica12_Character_064[] = { 12,  0,  0,  0,  0,  0,  0, 31,  0, 32,  0, 77,128, 83, 64, 81, 32, 81, 32, 73, 32, 38,160, 48, 64, 15,128,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_065[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 65,  0, 65,  0, 65,  0, 62,  0, 34,  0, 34,  0, 20,  0, 20,  0,  8,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_066[] = {  8,  0,  0,  0,  0,124, 66, 66, 66,124, 66, 66, 66,124,  0,  0,  0};
static const GLubyte Helvetica12_Character_067[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 30,  0, 33,  0, 64,  0, 64,  0, 64,  0, 64,  0, 64,  0, 33,  0, 30,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_068[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0,124,  0, 66,  0, 65,  0, 65,  0, 65,  0, 65,  0, 65,  0, 66,  0,124,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_069[] = {  8,  0,  0,  0,  0,126, 64, 64, 64,126, 64, 64, 64,126,  0,  0,  0};
static const GLubyte Helvetica12_Character_070[] = {  8,  0,  0,  0,  0, 64, 64, 64, 64,124, 64, 64, 64,126,  0,  0,  0};
static const GLubyte Helvetica12_Character_071[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 29,  0, 35,  0, 65,  0, 65,  0, 71,  0, 64,  0, 64,  0, 33,  0, 30,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_072[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 65,  0, 65,  0, 65,  0, 65,  0,127,  0, 65,  0, 65,  0, 65,  0, 65,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_073[] = {  3,  0,  0,  0,  0, 64, 64, 64, 64, 64, 64, 64, 64, 64,  0,  0,  0};
static const GLubyte Helvetica12_Character_074[] = {  7,  0,  0,  0,  0, 56, 68, 68,  4,  4,  4,  4,  4,  4,  0,  0,  0};
static const GLubyte Helvetica12_Character_075[] = {  8,  0,  0,  0,  0, 65, 66, 68, 72,112, 80, 72, 68, 66,  0,  0,  0};
static const GLubyte Helvetica12_Character_076[] = {  7,  0,  0,  0,  0,124, 64, 64, 64, 64, 64, 64, 64, 64,  0,  0,  0};
static const GLubyte Helvetica12_Character_077[] = { 11,  0,  0,  0,  0,  0,  0,  0,  0, 68, 64, 68, 64, 74, 64, 74, 64, 81, 64, 81, 64, 96,192, 96,192, 64, 64,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_078[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 65,  0, 67,  0, 69,  0, 69,  0, 73,  0, 81,  0, 81,  0, 97,  0, 65,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_079[] = { 10,  0,  0,  0,  0,  0,  0,  0,  0, 30,  0, 33,  0, 64,128, 64,128, 64,128, 64,128, 64,128, 33,  0, 30,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_080[] = {  8,  0,  0,  0,  0, 64, 64, 64, 64,124, 66, 66, 66,124,  0,  0,  0};
static const GLubyte Helvetica12_Character_081[] = { 10,  0,  0,  0,  0,  0,  0,  0,  0, 30,128, 33,  0, 66,128, 68,128, 64,128, 64,128, 64,128, 33,  0, 30,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_082[] = {  8,  0,  0,  0,  0, 66, 66, 66, 68,124, 66, 66, 66,124,  0,  0,  0};
static const GLubyte Helvetica12_Character_083[] = {  8,  0,  0,  0,  0, 60, 66, 66,  2, 12, 48, 64, 66, 60,  0,  0,  0};
static const GLubyte Helvetica12_Character_084[] = {  7,  0,  0,  0,  0, 16, 16, 16, 16, 16, 16, 16, 16,254,  0,  0,  0};
static const GLubyte Helvetica12_Character_085[] = {  8,  0,  0,  0,  0, 60, 66, 66, 66, 66, 66, 66, 66, 66,  0,  0,  0};
static const GLubyte Helvetica12_Character_086[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0,  8,  0,  8,  0, 20,  0, 20,  0, 34,  0, 34,  0, 34,  0, 65,  0, 65,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_087[] = { 11,  0,  0,  0,  0,  0,  0,  0,  0, 17,  0, 17,  0, 17,  0, 42,128, 42,128, 36,128, 68, 64, 68, 64, 68, 64,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_088[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 65,  0, 34,  0, 34,  0, 20,  0,  8,  0, 20,  0, 34,  0, 34,  0, 65,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_089[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0,  8,  0,  8,  0,  8,  0,  8,  0, 20,  0, 34,  0, 34,  0, 65,  0, 65,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_090[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0,127,  0, 64,  0, 32,  0, 16,  0,  8,  0,  4,  0,  2,  0,  1,  0,127,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_091[] = {  3,  0, 96, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 96,  0,  0,  0};
static const GLubyte Helvetica12_Character_092[] = {  4,  0,  0,  0,  0, 16, 16, 32, 32, 32, 64, 64,128,128,  0,  0,  0};
static const GLubyte Helvetica12_Character_093[] = {  3,  0,192, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,192,  0,  0,  0};
static const GLubyte Helvetica12_Character_094[] = {  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,136, 80, 32,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_095[] = {  7,  0,  0,254,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_096[] = {  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,192,128, 64,  0,  0,  0};
static const GLubyte Helvetica12_Character_097[] = {  7,  0,  0,  0,  0, 58, 68, 68, 60,  4, 68, 56,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_098[] = {  7,  0,  0,  0,  0, 88,100, 68, 68, 68,100, 88, 64, 64,  0,  0,  0};
static const GLubyte Helvetica12_Character_099[] = {  7,  0,  0,  0,  0, 56, 68, 64, 64, 64, 68, 56,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_100[] = {  7,  0,  0,  0,  0, 52, 76, 68, 68, 68, 76, 52,  4,  4,  0,  0,  0};
static const GLubyte Helvetica12_Character_101[] = {  7,  0,  0,  0,  0, 56, 68, 64,124, 68, 68, 56,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_102[] = {  3,  0,  0,  0,  0, 64, 64, 64, 64, 64, 64,224, 64, 48,  0,  0,  0};
static const GLubyte Helvetica12_Character_103[] = {  7,  0, 56, 68,  4, 52, 76, 68, 68, 68, 76, 52,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_104[] = {  7,  0,  0,  0,  0, 68, 68, 68, 68, 68,100, 88, 64, 64,  0,  0,  0};
static const GLubyte Helvetica12_Character_105[] = {  3,  0,  0,  0,  0, 64, 64, 64, 64, 64, 64, 64,  0, 64,  0,  0,  0};
static const GLubyte Helvetica12_Character_106[] = {  3,  0,128, 64, 64, 64, 64, 64, 64, 64, 64, 64,  0, 64,  0,  0,  0};
static const GLubyte Helvetica12_Character_107[] = {  6,  0,  0,  0,  0, 68, 72, 80, 96, 96, 80, 72, 64, 64,  0,  0,  0};
static const GLubyte Helvetica12_Character_108[] = {  3,  0,  0,  0,  0, 64, 64, 64, 64, 64, 64, 64, 64, 64,  0,  0,  0};
static const GLubyte Helvetica12_Character_109[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 73,  0, 73,  0, 73,  0, 73,  0, 73,  0,109,  0, 82,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_110[] = {  7,  0,  0,  0,  0, 68, 68, 68, 68, 68,100, 88,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_111[] = {  7,  0,  0,  0,  0, 56, 68, 68, 68, 68, 68, 56,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_112[] = {  7,  0, 64, 64, 64, 88,100, 68, 68, 68,100, 88,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_113[] = {  7,  0,  4,  4,  4, 52, 76, 68, 68, 68, 76, 52,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_114[] = {  4,  0,  0,  0,  0, 64, 64, 64, 64, 64, 96, 80,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_115[] = {  6,  0,  0,  0,  0, 48, 72,  8, 48, 64, 72, 48,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_116[] = {  3,  0,  0,  0,  0, 96, 64, 64, 64, 64, 64,224, 64, 64,  0,  0,  0};
static const GLubyte Helvetica12_Character_117[] = {  7,  0,  0,  0,  0, 52, 76, 68, 68, 68, 68, 68,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_118[] = {  7,  0,  0,  0,  0, 16, 16, 40, 40, 68, 68, 68,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_119[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 34,  0, 34,  0, 85,  0, 73,  0, 73,  0,136,128,136,128,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_120[] = {  6,  0,  0,  0,  0,132,132, 72, 48, 48, 72,132,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_121[] = {  7,  0, 64, 32, 16, 16, 40, 40, 72, 68, 68, 68,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_122[] = {  6,  0,  0,  0,  0,120, 64, 32, 32, 16,  8,120,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_123[] = {  4,  0, 48, 64, 64, 64, 64, 64,128, 64, 64, 64, 64, 48,  0,  0,  0};
static const GLubyte Helvetica12_Character_124[] = {  3,  0, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,  0,  0,  0};
static const GLubyte Helvetica12_Character_125[] = {  4,  0,192, 32, 32, 32, 32, 32, 16, 32, 32, 32, 32,192,  0,  0,  0};
static const GLubyte Helvetica12_Character_126[] = {  7,  0,  0,  0,  0,  0,  0,  0,152,100,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_127[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_128[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_129[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_130[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_131[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_132[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_133[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_134[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_135[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_136[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_137[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_138[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_139[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_140[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_141[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_142[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_143[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_144[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_145[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_146[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_147[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_148[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_149[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_150[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_151[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_152[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_153[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_154[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_155[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_156[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_157[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_158[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_159[] = { 12,  0,  0,  0,  0,  0,  0,  0,  0, 85,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 65,  0,  0,  0, 85,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_160[] = {  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_161[] = {  3,  0, 64, 64, 64, 64, 64, 64, 64, 64,  0, 64,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_162[] = {  7,  0,  0,  0, 32, 56,100, 80, 80, 80, 84, 56,  8,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_163[] = {  7,  0,  0,  0,  0, 88, 36, 16, 16,120, 32, 32, 36, 24,  0,  0,  0};
static const GLubyte Helvetica12_Character_164[] = {  7,  0,  0,  0,  0,  0,132,120, 72, 72,120,132,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_165[] = {  7,  0,  0,  0,  0, 16, 16,124, 16,124, 16, 40, 68, 68,  0,  0,  0};
static const GLubyte Helvetica12_Character_166[] = {  3,  0,  0, 64, 64, 64, 64,  0,  0,  0, 64, 64, 64, 64,  0,  0,  0};
static const GLubyte Helvetica12_Character_167[] = {  6,  0,112,136,  8, 48, 72,136,136,144, 96,128,136,112,  0,  0,  0};
static const GLubyte Helvetica12_Character_168[] = {  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,160,  0,  0,  0};
static const GLubyte Helvetica12_Character_169[] = { 11,  0,  0,  0,  0,  0,  0,  0,  0, 31,  0, 32,128, 78, 64, 81, 64, 80, 64, 81, 64, 78, 64, 32,128, 31,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_170[] = {  5,  0,  0,  0,  0,  0,  0,  0,  0,112,  0, 80, 16,112,  0,  0,  0};
static const GLubyte Helvetica12_Character_171[] = {  7,  0,  0,  0,  0,  0, 20, 40, 80, 40, 20,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_172[] = {  8,  0,  0,  0,  0,  0,  0,  2,  2,  2,126,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_173[] = {  5,  0,  0,  0,  0,  0,  0,  0,240,  0,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_174[] = { 11,  0,  0,  0,  0,  0,  0,  0,  0, 31,  0, 32,128, 74, 64, 74, 64, 76, 64, 74, 64, 78, 64, 32,128, 31,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_175[] = {  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,240,  0,  0,  0};
static const GLubyte Helvetica12_Character_176[] = {  5,  0,  0,  0,  0,  0,  0,  0,  0, 96,144,144, 96,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_177[] = {  7,  0,  0,  0,  0,124,  0, 16, 16,124, 16, 16,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_178[] = {  4,  0,  0,  0,  0,  0,  0,  0,240, 64, 32,144, 96,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_179[] = {  4,  0,  0,  0,  0,  0,  0,  0,192, 32, 64, 32,224,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_180[] = {  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,128, 64,  0,  0};
static const GLubyte Helvetica12_Character_181[] = {  7,  0, 64, 64, 64,116, 76, 68, 68, 68, 68, 68,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_182[] = {  7,  0, 40, 40, 40, 40, 40, 40,104,232,232,232,104, 60,  0,  0,  0};
static const GLubyte Helvetica12_Character_183[] = {  3,  0,  0,  0,  0,  0,  0,  0, 64,  0,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_184[] = {  3,  0,192, 32, 32, 64,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_185[] = {  4,  0,  0,  0,  0,  0,  0,  0, 32, 32, 32, 96, 32,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_186[] = {  5,  0,  0,  0,  0,  0,  0,  0,  0,112,  0,112, 80,112,  0,  0,  0};
static const GLubyte Helvetica12_Character_187[] = {  7,  0,  0,  0,  0,  0, 80, 40, 20, 40, 80,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_188[] = { 10,  0,  0,  0,  0,  0,  0,  0,  0, 65,  0, 39,128, 21,  0, 19,  0, 73,  0, 68,  0, 68,  0,194,  0, 65,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_189[] = { 10,  0,  0,  0,  0,  0,  0,  0,  0, 71,128, 34,  0, 17,  0, 20,128, 75,  0, 72,  0, 68,  0,194,  0, 65,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_190[] = { 10,  0,  0,  0,  0,  0,  0,  0,  0, 33,  0, 23,128, 21,  0, 11,  0,201,  0, 36,  0, 68,  0, 34,  0,225,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_191[] = {  7,  0, 56, 68, 68, 32, 32, 16, 16,  0, 16,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_192[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 65,  0, 65,  0, 65,  0, 62,  0, 34,  0, 34,  0, 20,  0,  8,  0,  8,  0,  0,  0,  8,  0, 16,  0};
static const GLubyte Helvetica12_Character_193[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 65,  0, 65,  0, 65,  0, 62,  0, 34,  0, 34,  0, 20,  0,  8,  0,  8,  0,  0,  0,  8,  0,  4,  0};
static const GLubyte Helvetica12_Character_194[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 65,  0, 65,  0, 65,  0, 62,  0, 34,  0, 34,  0, 20,  0,  8,  0,  8,  0,  0,  0, 20,  0,  8,  0};
static const GLubyte Helvetica12_Character_195[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 65,  0, 65,  0, 65,  0, 62,  0, 34,  0, 34,  0, 20,  0,  8,  0,  8,  0,  0,  0, 20,  0, 10,  0};
static const GLubyte Helvetica12_Character_196[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 65,  0, 65,  0, 65,  0, 62,  0, 34,  0, 34,  0, 20,  0,  8,  0,  8,  0,  0,  0, 20,  0,  0,  0};
static const GLubyte Helvetica12_Character_197[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 65,  0, 65,  0, 65,  0, 62,  0, 34,  0, 34,  0, 20,  0,  8,  0,  8,  0,  8,  0, 20,  0,  8,  0};
static const GLubyte Helvetica12_Character_198[] = { 11,  0,  0,  0,  0,  0,  0,  0,  0, 71,192, 68,  0, 68,  0, 60,  0, 39,192, 36,  0, 20,  0, 20,  0, 15,192,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_199[] = {  9,  0,  0, 24,  0,  4,  0,  4,  0, 30,  0, 33,  0, 64,  0, 64,  0, 64,  0, 64,  0, 64,  0, 33,  0, 30,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_200[] = {  8,  0,  0,  0,  0,126, 64, 64, 64,126, 64, 64, 64,126,  0,  8, 16};
static const GLubyte Helvetica12_Character_201[] = {  8,  0,  0,  0,  0,126, 64, 64, 64,126, 64, 64, 64,126,  0,  8,  4};
static const GLubyte Helvetica12_Character_202[] = {  8,  0,  0,  0,  0,126, 64, 64, 64,126, 64, 64, 64,126,  0, 20,  8};
static const GLubyte Helvetica12_Character_203[] = {  8,  0,  0,  0,  0,126, 64, 64, 64,126, 64, 64, 64,126,  0, 20,  0};
static const GLubyte Helvetica12_Character_204[] = {  3,  0,  0,  0,  0, 64, 64, 64, 64, 64, 64, 64, 64, 64,  0, 64,128};
static const GLubyte Helvetica12_Character_205[] = {  3,  0,  0,  0,  0, 64, 64, 64, 64, 64, 64, 64, 64, 64,  0, 64, 32};
static const GLubyte Helvetica12_Character_206[] = {  3,  0,  0,  0,  0, 64, 64, 64, 64, 64, 64, 64, 64, 64,  0,160, 64};
static const GLubyte Helvetica12_Character_207[] = {  3,  0,  0,  0,  0, 64, 64, 64, 64, 64, 64, 64, 64, 64,  0,160,  0};
static const GLubyte Helvetica12_Character_208[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0,124,  0, 66,  0, 65,  0, 65,  0,241,  0, 65,  0, 65,  0, 66,  0,124,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_209[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0, 65,  0, 67,  0, 69,  0, 69,  0, 73,  0, 81,  0, 81,  0, 97,  0, 65,  0,  0,  0, 20,  0, 10,  0};
static const GLubyte Helvetica12_Character_210[] = { 10,  0,  0,  0,  0,  0,  0,  0,  0, 30,  0, 33,  0, 64,128, 64,128, 64,128, 64,128, 64,128, 33,  0, 30,  0,  0,  0,  4,  0,  8,  0};
static const GLubyte Helvetica12_Character_211[] = { 10,  0,  0,  0,  0,  0,  0,  0,  0, 30,  0, 33,  0, 64,128, 64,128, 64,128, 64,128, 64,128, 33,  0, 30,  0,  0,  0,  4,  0,  2,  0};
static const GLubyte Helvetica12_Character_212[] = { 10,  0,  0,  0,  0,  0,  0,  0,  0, 30,  0, 33,  0, 64,128, 64,128, 64,128, 64,128, 64,128, 33,  0, 30,  0,  0,  0, 10,  0,  4,  0};
static const GLubyte Helvetica12_Character_213[] = { 10,  0,  0,  0,  0,  0,  0,  0,  0, 30,  0, 33,  0, 64,128, 64,128, 64,128, 64,128, 64,128, 33,  0, 30,  0,  0,  0, 20,  0, 10,  0};
static const GLubyte Helvetica12_Character_214[] = { 10,  0,  0,  0,  0,  0,  0,  0,  0, 30,  0, 33,  0, 64,128, 64,128, 64,128, 64,128, 64,128, 33,  0, 30,  0,  0,  0, 18,  0,  0,  0};
static const GLubyte Helvetica12_Character_215[] = {  7,  0,  0,  0,  0,  0, 68, 40, 16, 40, 68,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_216[] = { 10,  0,  0,  0,  0,  0,  0,128,  0, 94,  0, 33,  0, 80,128, 72,128, 68,128, 68,128, 66,128, 33,  0, 30,128,  0, 64,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_217[] = {  8,  0,  0,  0,  0, 60, 66, 66, 66, 66, 66, 66, 66, 66,  0,  8, 16};
static const GLubyte Helvetica12_Character_218[] = {  8,  0,  0,  0,  0, 60, 66, 66, 66, 66, 66, 66, 66, 66,  0,  8,  4};
static const GLubyte Helvetica12_Character_219[] = {  8,  0,  0,  0,  0, 60, 66, 66, 66, 66, 66, 66, 66, 66,  0, 20,  8};
static const GLubyte Helvetica12_Character_220[] = {  8,  0,  0,  0,  0, 60, 66, 66, 66, 66, 66, 66, 66, 66,  0, 36,  0};
static const GLubyte Helvetica12_Character_221[] = {  9,  0,  0,  0,  0,  0,  0,  0,  0,  8,  0,  8,  0,  8,  0,  8,  0, 20,  0, 34,  0, 34,  0, 65,  0, 65,  0,  0,  0,  8,  0,  4,  0};
static const GLubyte Helvetica12_Character_222[] = {  8,  0,  0,  0,  0, 64, 64,124, 66, 66, 66,124, 64, 64,  0,  0,  0};
static const GLubyte Helvetica12_Character_223[] = {  7,  0,  0,  0,  0, 88, 68, 68, 68, 88, 68, 68, 68, 56,  0,  0,  0};
static const GLubyte Helvetica12_Character_224[] = {  7,  0,  0,  0,  0, 58, 68, 68, 60,  4, 68, 56,  0,  8, 16,  0,  0};
static const GLubyte Helvetica12_Character_225[] = {  7,  0,  0,  0,  0, 58, 68, 68, 60,  4, 68, 56,  0, 16,  8,  0,  0};
static const GLubyte Helvetica12_Character_226[] = {  7,  0,  0,  0,  0, 58, 68, 68, 60,  4, 68, 56,  0, 40, 16,  0,  0};
static const GLubyte Helvetica12_Character_227[] = {  7,  0,  0,  0,  0, 58, 68, 68, 60,  4, 68, 56,  0, 40, 20,  0,  0};
static const GLubyte Helvetica12_Character_228[] = {  7,  0,  0,  0,  0, 58, 68, 68, 60,  4, 68, 56,  0, 40,  0,  0,  0};
static const GLubyte Helvetica12_Character_229[] = {  7,  0,  0,  0,  0, 58, 68, 68, 60,  4, 68, 56, 24, 36, 24,  0,  0};
static const GLubyte Helvetica12_Character_230[] = { 11,  0,  0,  0,  0,  0,  0,  0,  0, 59,128, 68, 64, 68,  0, 63,192,  4, 64, 68, 64, 59,128,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_231[] = {  7,  0, 48,  8, 16, 56, 68, 64, 64, 64, 68, 56,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_232[] = {  7,  0,  0,  0,  0, 56, 68, 64,124, 68, 68, 56,  0, 16, 32,  0,  0};
static const GLubyte Helvetica12_Character_233[] = {  7,  0,  0,  0,  0, 56, 68, 64,124, 68, 68, 56,  0, 16,  8,  0,  0};
static const GLubyte Helvetica12_Character_234[] = {  7,  0,  0,  0,  0, 56, 68, 64,124, 68, 68, 56,  0, 40, 16,  0,  0};
static const GLubyte Helvetica12_Character_235[] = {  7,  0,  0,  0,  0, 56, 68, 64,124, 68, 68, 56,  0, 40,  0,  0,  0};
static const GLubyte Helvetica12_Character_236[] = {  3,  0,  0,  0,  0, 64, 64, 64, 64, 64, 64, 64,  0, 64,128,  0,  0};
static const GLubyte Helvetica12_Character_237[] = {  3,  0,  0,  0,  0, 64, 64, 64, 64, 64, 64, 64,  0, 64, 32,  0,  0};
static const GLubyte Helvetica12_Character_238[] = {  3,  0,  0,  0,  0, 64, 64, 64, 64, 64, 64, 64,  0,160, 64,  0,  0};
static const GLubyte Helvetica12_Character_239[] = {  3,  0,  0,  0,  0, 64, 64, 64, 64, 64, 64, 64,  0,160,  0,  0,  0};
static const GLubyte Helvetica12_Character_240[] = {  7,  0,  0,  0,  0, 56, 68, 68, 68, 68, 60,  4, 40, 24, 52,  0,  0};
static const GLubyte Helvetica12_Character_241[] = {  7,  0,  0,  0,  0, 68, 68, 68, 68, 68,100, 88,  0, 40, 20,  0,  0};
static const GLubyte Helvetica12_Character_242[] = {  7,  0,  0,  0,  0, 56, 68, 68, 68, 68, 68, 56,  0, 16, 32,  0,  0};
static const GLubyte Helvetica12_Character_243[] = {  7,  0,  0,  0,  0, 56, 68, 68, 68, 68, 68, 56,  0, 16,  8,  0,  0};
static const GLubyte Helvetica12_Character_244[] = {  7,  0,  0,  0,  0, 56, 68, 68, 68, 68, 68, 56,  0, 40, 16,  0,  0};
static const GLubyte Helvetica12_Character_245[] = {  7,  0,  0,  0,  0, 56, 68, 68, 68, 68, 68, 56,  0, 40, 20,  0,  0};
static const GLubyte Helvetica12_Character_246[] = {  7,  0,  0,  0,  0, 56, 68, 68, 68, 68, 68, 56,  0, 40,  0,  0,  0};
static const GLubyte Helvetica12_Character_247[] = {  7,  0,  0,  0,  0,  0, 16,  0,124,  0, 16,  0,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_248[] = {  7,  0,  0,  0,  0,184, 68,100, 84, 76, 68, 58,  0,  0,  0,  0,  0};
static const GLubyte Helvetica12_Character_249[] = {  7,  0,  0,  0,  0, 52, 76, 68, 68, 68, 68, 68,  0, 16, 32,  0,  0};
static const GLubyte Helvetica12_Character_250[] = {  7,  0,  0,  0,  0, 52, 76, 68, 68, 68, 68, 68,  0, 16,  8,  0,  0};
static const GLubyte Helvetica12_Character_251[] = {  7,  0,  0,  0,  0, 52, 76, 68, 68, 68, 68, 68,  0, 40, 16,  0,  0};
static const GLubyte Helvetica12_Character_252[] = {  7,  0,  0,  0,  0, 52, 76, 68, 68, 68, 68, 68,  0, 40,  0,  0,  0};
static const GLubyte Helvetica12_Character_253[] = {  7,  0, 64, 32, 16, 16, 40, 40, 72, 68, 68, 68,  0, 16,  8,  0,  0};
static const GLubyte Helvetica12_Character_254[] = {  7,  0, 64, 64, 64, 88,100, 68, 68, 68,100, 88, 64, 64,  0,  0,  0};
static const GLubyte Helvetica12_Character_255[] = {  7,  0, 96, 16, 16, 16, 24, 40, 40, 36, 68, 68,  0, 40,  0,  0,  0};

/* The font characters mapping: */
static const GLubyte* Helvetica12_Character_Map[] = {Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,
                                                     Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,
                                                     Helvetica12_Character_032,Helvetica12_Character_033,Helvetica12_Character_034,Helvetica12_Character_035,Helvetica12_Character_036,Helvetica12_Character_037,Helvetica12_Character_038,Helvetica12_Character_039,Helvetica12_Character_040,Helvetica12_Character_041,Helvetica12_Character_042,Helvetica12_Character_043,Helvetica12_Character_044,Helvetica12_Character_045,Helvetica12_Character_046,Helvetica12_Character_047,
                                                     Helvetica12_Character_048,Helvetica12_Character_049,Helvetica12_Character_050,Helvetica12_Character_051,Helvetica12_Character_052,Helvetica12_Character_053,Helvetica12_Character_054,Helvetica12_Character_055,Helvetica12_Character_056,Helvetica12_Character_057,Helvetica12_Character_058,Helvetica12_Character_059,Helvetica12_Character_060,Helvetica12_Character_061,Helvetica12_Character_062,Helvetica12_Character_063,
                                                     Helvetica12_Character_064,Helvetica12_Character_065,Helvetica12_Character_066,Helvetica12_Character_067,Helvetica12_Character_068,Helvetica12_Character_069,Helvetica12_Character_070,Helvetica12_Character_071,Helvetica12_Character_072,Helvetica12_Character_073,Helvetica12_Character_074,Helvetica12_Character_075,Helvetica12_Character_076,Helvetica12_Character_077,Helvetica12_Character_078,Helvetica12_Character_079,
                                                     Helvetica12_Character_080,Helvetica12_Character_081,Helvetica12_Character_082,Helvetica12_Character_083,Helvetica12_Character_084,Helvetica12_Character_085,Helvetica12_Character_086,Helvetica12_Character_087,Helvetica12_Character_088,Helvetica12_Character_089,Helvetica12_Character_090,Helvetica12_Character_091,Helvetica12_Character_092,Helvetica12_Character_093,Helvetica12_Character_094,Helvetica12_Character_095,
                                                     Helvetica12_Character_096,Helvetica12_Character_097,Helvetica12_Character_098,Helvetica12_Character_099,Helvetica12_Character_100,Helvetica12_Character_101,Helvetica12_Character_102,Helvetica12_Character_103,Helvetica12_Character_104,Helvetica12_Character_105,Helvetica12_Character_106,Helvetica12_Character_107,Helvetica12_Character_108,Helvetica12_Character_109,Helvetica12_Character_110,Helvetica12_Character_111,
                                                     Helvetica12_Character_112,Helvetica12_Character_113,Helvetica12_Character_114,Helvetica12_Character_115,Helvetica12_Character_116,Helvetica12_Character_117,Helvetica12_Character_118,Helvetica12_Character_119,Helvetica12_Character_120,Helvetica12_Character_121,Helvetica12_Character_122,Helvetica12_Character_123,Helvetica12_Character_124,Helvetica12_Character_125,Helvetica12_Character_126,Helvetica12_Character_032,
                                                     Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,
                                                     Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,Helvetica12_Character_032,
                                                     Helvetica12_Character_160,Helvetica12_Character_161,Helvetica12_Character_162,Helvetica12_Character_163,Helvetica12_Character_164,Helvetica12_Character_165,Helvetica12_Character_166,Helvetica12_Character_167,Helvetica12_Character_168,Helvetica12_Character_169,Helvetica12_Character_170,Helvetica12_Character_171,Helvetica12_Character_172,Helvetica12_Character_173,Helvetica12_Character_174,Helvetica12_Character_175,
                                                     Helvetica12_Character_176,Helvetica12_Character_177,Helvetica12_Character_178,Helvetica12_Character_179,Helvetica12_Character_180,Helvetica12_Character_181,Helvetica12_Character_182,Helvetica12_Character_183,Helvetica12_Character_184,Helvetica12_Character_185,Helvetica12_Character_186,Helvetica12_Character_187,Helvetica12_Character_188,Helvetica12_Character_189,Helvetica12_Character_190,Helvetica12_Character_191,
                                                     Helvetica12_Character_192,Helvetica12_Character_193,Helvetica12_Character_194,Helvetica12_Character_195,Helvetica12_Character_196,Helvetica12_Character_197,Helvetica12_Character_198,Helvetica12_Character_199,Helvetica12_Character_200,Helvetica12_Character_201,Helvetica12_Character_202,Helvetica12_Character_203,Helvetica12_Character_204,Helvetica12_Character_205,Helvetica12_Character_206,Helvetica12_Character_207,
                                                     Helvetica12_Character_208,Helvetica12_Character_209,Helvetica12_Character_210,Helvetica12_Character_211,Helvetica12_Character_212,Helvetica12_Character_213,Helvetica12_Character_214,Helvetica12_Character_215,Helvetica12_Character_216,Helvetica12_Character_217,Helvetica12_Character_218,Helvetica12_Character_219,Helvetica12_Character_220,Helvetica12_Character_221,Helvetica12_Character_222,Helvetica12_Character_223,
                                                     Helvetica12_Character_224,Helvetica12_Character_225,Helvetica12_Character_226,Helvetica12_Character_227,Helvetica12_Character_228,Helvetica12_Character_229,Helvetica12_Character_230,Helvetica12_Character_231,Helvetica12_Character_232,Helvetica12_Character_233,Helvetica12_Character_234,Helvetica12_Character_235,Helvetica12_Character_236,Helvetica12_Character_237,Helvetica12_Character_238,Helvetica12_Character_239,
                                                     Helvetica12_Character_240,Helvetica12_Character_241,Helvetica12_Character_242,Helvetica12_Character_243,Helvetica12_Character_244,Helvetica12_Character_245,Helvetica12_Character_246,Helvetica12_Character_247,Helvetica12_Character_248,Helvetica12_Character_249,Helvetica12_Character_250,Helvetica12_Character_251,Helvetica12_Character_252,Helvetica12_Character_253,Helvetica12_Character_254,Helvetica12_Character_255,NULL};

/* The font structure: */
static const NBody_SFG_Font nbody_fgFontHelvetica12 = { "-adobe-helvetica-medium-r-normal--12-120-75-75-p-67-iso8859-1", 256, 16, Helvetica12_Character_Map, 0, 4 };

void nbody_glutBitmapStringHelvetica(const unsigned char* string)
{
    unsigned char c;
    float x = 0.0f;
    const NBody_SFG_Font* font = &nbody_fgFontHelvetica12;

    glPushClientAttrib(GL_CLIENT_PIXEL_STORE_BIT);
    glPixelStorei(GL_UNPACK_SWAP_BYTES,  GL_FALSE);
    glPixelStorei(GL_UNPACK_LSB_FIRST,   GL_FALSE);
    glPixelStorei(GL_UNPACK_ROW_LENGTH,  0       );
    glPixelStorei(GL_UNPACK_SKIP_ROWS,   0       );
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0       );
    glPixelStorei(GL_UNPACK_ALIGNMENT,   1       );

    /*
     * Step through the string, drawing each character.
     * A newline will simply translate the next character's insertion
     * point back to the start of the line and down one line.
     */
    while ((c = *string++))
    {
        if (c == '\n')
        {
            glBitmap (0, 0, 0, 0, -x, (float) -font->Height, NULL);
            x = 0.0f;
        }
        else  /* Not an EOL, draw the bitmap character */
        {
            const GLubyte* face = font->Characters[c];

            glBitmap(
                face[0], font->Height,     /* Bitmap's width and height    */
                font->xorig, font->yorig,  /* The origin in the font glyph */
                (float) face[0], 0.0,      /* The raster advance; inc. x,y */
                (face + 1)                 /* The packed bitmap data...    */
            );

            x += (float)face[0];
        }
    }

    glPopClientAttrib();
}


