/*
Copyright 2008-2010 Travis Desell, Dave Przybylo, Nathan Cole, Matthew
Arsenault, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <windows.h>
#include <stdio.h>

#include <boinc/str_util.h>

extern "C" int main(int argc, const char* argv[]);

/* FIXME: undefined reference to _set_invalid_parameter_handler */
void AppInvalidParameterHandler(const wchar_t* expression,
                                const wchar_t* function,
                                const wchar_t* file,
                                unsigned int line,
                                uintptr_t pReserved )
{
    fprintf(stderr, "Invalid parameter detected in function %ls. File: %ls Line: %d\n", function, file, line);
    fprintf(stderr, "Expression: %ls\n", expression);
    // Cause a Debug Breakpoint.
    DebugBreak();
}

int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR Args, int WinMode)
{
    LPSTR command_line;
    char* argv[100];
    int argc;

    //_set_invalid_parameter_handler(AppInvalidParameterHandler);

    command_line = GetCommandLine();
    argc = parse_command_line( command_line, argv );
    return main(argc, (const char**) argv);
}


