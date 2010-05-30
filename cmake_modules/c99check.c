/*
 * Multi precision toolbox for Scilab
 * Copyright (C) 2009 - Jonathan Blanchard
 *
 * This file must be used under the terms of the CeCILL.
 * This source file is licensed as described in the file COPYING, which
 * you should have received as part of this distribution.  The terms
 * are also available at
 * http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
 */

/* Source that check for C99 compiler support. */

#if __STDC_VERSION__ >= 199901L
    #error "No C99 support detected."
#endif


int main(int argc, char *argv[])
{
    /* Dummy stuff just in case. */
    if( argc > 100 )
        return 1;

    return 0;
}
