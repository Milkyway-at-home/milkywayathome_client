/****************************************************************************
 *                                                                          *
 *  Copyright (C) 2010 Shane Reilly and Rensselaer Polytechnic Institute.   *
 *                                                                          *
 *  This file is part of the Light Modeling Library (LModL).                *
 *                                                                          *
 *  This library is free software: you can redistribute it and/or modify    *
 *  it under the terms of the GNU General Public License as published by    *
 *  the Free Software Foundation, either version 3 of the License, or       *
 *  (at your option) any later version.                                     *
 *                                                                          *
 *  This library is distributed in the hope that it will be useful,         *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *  GNU General Public License for more details.                            *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License       *
 *  along with this library. If not, see <http://www.gnu.org/licenses/>.    *
 *                                                                          *
 *  Shane Reilly                                                            *
 *  reills2@cs.rpi.edu                                                      *
 *                                                                          *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

//#define TEST_MODE

int isBigEndian()
{
    int endianTest = 1;
    if ( *( (char*) &endianTest ) == 1 )
        return 0;
    else
        return 1;
}

inline void endianSwapInt16( unsigned long i, char* out )
{

    unsigned char* ia = (unsigned char*) &i;

    out[0] = ia[1];
    out[1] = ia[0];

}

inline void endianSwapInt32( unsigned long i, char* out )
{

    unsigned char* ia = (unsigned char*) &i;

    out[0] = ia[3];
    out[1] = ia[2];
    out[2] = ia[1];
    out[3] = ia[0];

}

inline void endianSwapInt64( unsigned long i, char* out )
{

    unsigned char* ia = (unsigned char*) &i;

    out[0] = ia[7];
    out[1] = ia[6];
    out[2] = ia[5];
    out[3] = ia[4];
    out[4] = ia[3];
    out[5] = ia[2];
    out[6] = ia[1];
    out[7] = ia[0];

}

inline void endianSwapfloat( float d, char* out )
{

    unsigned char* da = (unsigned char*) &d;

    out[0] = da[7];
    out[1] = da[6];
    out[2] = da[5];
    out[3] = da[4];
    out[4] = da[3];
    out[5] = da[2];
    out[6] = da[1];
    out[7] = da[0];

}

inline void getBytesInt16( unsigned long i, char* out )
{

    unsigned char* ia = (unsigned char*) &i;

    out[0] = ia[0];
    out[1] = ia[1];

}

inline void getBytesInt32( unsigned long i, char* out )
{

    unsigned char* ia = (unsigned char*) &i;

    out[0] = ia[0];
    out[1] = ia[1];
    out[2] = ia[2];
    out[3] = ia[3];

}

inline void getBytesInt64( unsigned long i, char* out )
{

    unsigned char* ia = (unsigned char*) &i;

    out[0] = ia[0];
    out[1] = ia[1];
    out[2] = ia[2];
    out[3] = ia[3];
    out[4] = ia[4];
    out[5] = ia[5];
    out[6] = ia[6];
    out[7] = ia[7];

}

inline void getBytesfloat( float d, char* out )
{

    unsigned char* da = (unsigned char*) &d;

    out[0] = da[0];
    out[1] = da[1];
    out[2] = da[2];
    out[3] = da[3];
    out[4] = da[4];
    out[5] = da[5];
    out[6] = da[6];
    out[7] = da[7];

}

inline void endianSwapFloat( float d, char* out )
{

    unsigned char* da = (unsigned char*) &d;

    out[0] = da[3];
    out[1] = da[2];
    out[2] = da[1];
    out[3] = da[0];

}

inline void getBytesFloat( float d, char* out )
{

    unsigned char* da = (unsigned char*) &d;

    out[0] = da[0];
    out[1] = da[1];
    out[2] = da[2];
    out[3] = da[3];

}

char getYN()
    // Returns 'y' or 'n' depending on which is pressed by user
{

    char c;
    do {

        c = tolower(getchar());
/*      if( !iscntrl(c) ) {
            putchar(c);
            putchar('\b');
        }
*/
    } while( c != 'y' && c != 'n' );
    printf("\n");

    return c;

}

inline int isNum( char c )
{

    switch( c ) {

    case '-':
    case '+':
    case '.':
    case 'E':
    case 'e':
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
        return 1;

    default:
        return 0;

    }

}

void usage()
{

    printf("Usage: ./text2bin input_file output_file -f|-d\n");
    printf("    -f  use floating point format (4 bytes)\n");
    printf("    -d  use float precision format (8 bytes)\n");
    exit(0);

}

int main( int args, char** argv )
{

#ifdef TEST_MODE
printf("DEBUG: checking arguments...\n");
#endif
    if( args<4 )
        usage();
    int useDblFlag;
    if( !strcmp(argv[3], "-f") )
        useDblFlag = 0;
    else if( !strcmp(argv[3], "-d") )
        useDblFlag = 1;
    else
        usage();
#ifdef TEST_MODE
printf("DEBUG: reading file...\n");
#endif
    FILE *inFile = fopen(argv[1], "r");
    if( inFile==NULL ) {
        fprintf(stderr, "Could not read '%s'.\n", argv[1]);
        exit(1);
    }
#ifdef TEST_MODE
printf("DEBUG: checking if write file exists...\n");
#endif
    // See if write-file exists
    FILE *outFile = fopen(argv[2], "r");
    if( outFile!=NULL ){
#ifdef TEST_MODE
printf("DEBUG: closing file...\n");
#endif
        fclose(outFile);
#ifdef TEST_MODE
printf("DEBUG: querying user...\n");
#endif
        printf("Overwrite file '%s'? (y/n): ", argv[2]);
        if( getYN()=='n' )
            return 0;
    }

#ifdef TEST_MODE
printf("DEBUG: opening write file...\n");
#endif
    outFile = fopen(argv[2], "wb");
    if( outFile==NULL ) {
        fprintf(stderr, "Could not open '%s' for writing.\n", argv[2]);
        exit(1);
    }

#ifdef TEST_MODE
printf("DEBUG: fetching infile character 1...\n");
#endif
    char c = fgetc(inFile);

    while( c!=EOF ) {

        // Check for number
#ifdef TEST_MODE
printf("DEBUG: see if input is a value\n");
#endif
        if( isNum(c) ) {

            // Store number string in 'value'
            char value[64];
            int i = 0;
            do {
                // Check for size
                if( i>=sizeof(value)-1 ) {
                    value[i++] = '\0';
                    fprintf(stderr, "Unable to parse data - missing delimiters: '%s'\n", value);
                    exit(1);
                }
                value[i++] = c;
                c = fgetc(inFile);
            } while( isNum(c) );
            value[i++] = '\0';

            // Determine if number is integer (4-bytes max in this version) or float type

#ifdef TEST_MODE
printf("DEBUG: check for integer/float type\n");
#endif
            char* iEndPtr;
            char* dEndPtr;

            unsigned long iConv = (unsigned long) strtol(value, &iEndPtr, 10);
            float dConv = strtod(value, &dEndPtr);

            if( iEndPtr==value+i-1 ) {
#ifdef TEST_MODE
printf("DEBUG: store integer value\n");
#endif
                /// TODO /// test this
                unsigned char da[8];
                // Value is an int
                if( isBigEndian() )
                    endianSwapInt32(dConv, da);
                else
                    getBytesInt32(dConv, da);
                // Write to file in binary format
                int i;
                for( i = 0; i<4; i++ )
                    if( fputc(da[i], outFile)==EOF ) {
                        fprintf(stderr, "Unable to write to file: '%s'\n", argv[2]);
                        exit(1);
                    }

            }
            else if( dEndPtr==value+i-1 ) {
                if( useDblFlag ) {
#ifdef TEST_MODE
printf("DEBUG: store float value\n");
#endif

                    char da[8];
                    // Value is a float
                    if( isBigEndian() )
                        endianSwapfloat(dConv, da);
                    else
                        getBytesfloat(dConv, da);
                    // Write to file in binary format
                    int i;
                    for( i = 0; i<8; i++ )
                        if( fputc(da[i], outFile)==EOF ) {
                            fprintf(stderr, "Unable to write to file: '%s'\n", argv[2]);
                            exit(1);
                        }

                }
                else {
#ifdef TEST_MODE
printf("DEBUG: store float value as a float\n");
#endif

                    char da[4];
                    // Value is a float
                    if( isBigEndian() )
                        endianSwapFloat((float) dConv, da);
                    else
                        getBytesFloat((float) dConv, da);
                    // Write to file in binary format
                    int i;
                    for( i = 0; i<4; i++ )
                        if( fputc(da[i], outFile)==EOF ) {
                            fprintf(stderr, "Unable to write to file: '%s'\n", argv[2]);
                            exit(1);
                        }

                }
            }
            else {
                fprintf(stderr, "Unable to parse data - value is not a valid int or float: '%s'\n", value);
                exit(1);
            }

        }

        // Make sure delimiters are plausible
#ifdef TEST_MODE
printf("DEBUG: check for eof\n");
#endif

        if( c==EOF )
            break;

#ifdef TEST_MODE
printf("DEBUG: test for delimiter\n");
#endif
        if( c!=' ' && !iscntrl(c) && !ispunct(c) ) {
            char str[64];
            str[0] = c;
            str[1] = '\0';
            fgets(str+1, sizeof(str)-1, inFile);
            fprintf(stderr, "Unable to parse data - string is not a valid int or float: '%s'\n", str);
            exit(1);
        }

#ifdef TEST_MODE
printf("DEBUG: input next character in file\n");
#endif

        c = fgetc(inFile);

    }

#ifdef TEST_MODE
printf("DEBUG: closing files\n");
#endif

    fclose(inFile);
    fclose(outFile);

    return 0;

}
