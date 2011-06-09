/*
Copyright (C) 2010  Matthew Arsenault

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


#include "milkyway_util.h"
#include "show_cal_types.h"

const char* showCALboolean(CALboolean x)
{
    switch (x)
    {
        case CAL_TRUE:
            return "CAL_TRUE";
        case CAL_FALSE:
            return "CAL_FALSE";
        default:
            warn("Trying to show unhandled CALbool: %d\n", x);
            return "Unknown CALbool";
    }
}

const char* showCALresult(CALresult x)
{
    switch (x)
    {
        case CAL_RESULT_OK:
            return "CAL_RESULT_OK";
        case CAL_RESULT_ERROR:
            return "CAL_RESULT_ERROR";
        case CAL_RESULT_INVALID_PARAMETER:
            return "CAL_RESULT_INVALID_PARAMETER";
        case CAL_RESULT_NOT_SUPPORTED:
            return "CAL_RESULT_NOT_SUPPORTED";
        case CAL_RESULT_ALREADY:
            return "CAL_RESULT_ALREADY";
        case CAL_RESULT_NOT_INITIALIZED:
            return "CAL_RESULT_NOT_INITIALIZED";
        case CAL_RESULT_BAD_HANDLE:
            return "CAL_RESULT_BAD_HANDLE";
        case CAL_RESULT_BAD_NAME_TYPE:
            return "CAL_RESULT_BAD_NAME_TYPE";
        case CAL_RESULT_PENDING:
            return "CAL_RESULT_PENDING";
        case CAL_RESULT_BUSY:
            return "CAL_RESULT_BUSY";
        case CAL_RESULT_WARNING:
            return "CAL_RESULT_WARNING";
        default:
            warn("Trying to show unhandled CALresult: %d\n", x);
            return "Unknown CALresult";
    }
}

const char* showCALtargetEnum(enum CALtargetEnum x)
{
    switch (x)
    {
        case CAL_TARGET_600:
            return "CAL_TARGET_600";
        case CAL_TARGET_610:
            return "CAL_TARGET_610";
        case CAL_TARGET_630:
            return "CAL_TARGET_630";
        case CAL_TARGET_670:
            return "CAL_TARGET_670";
        case CAL_TARGET_7XX:
            return "CAL_TARGET_7XX";
        case CAL_TARGET_770:
            return "CAL_TARGET_770";
        case CAL_TARGET_710:
            return "CAL_TARGET_710";
        case CAL_TARGET_730:
            return "CAL_TARGET_730";
        case CAL_TARGET_CYPRESS:
            return "CAL_TARGET_CYPRESS";
        case CAL_TARGET_JUNIPER:
            return "CAL_TARGET_JUNIPER";
        case CAL_TARGET_REDWOOD:
            return "CAL_TARGET_REDWOOD";
        case CAL_TARGET_CEDAR :
            return "CAL_TARGET_CEDAR";
        case CAL_TARGET_WRESTLER:
            return "CAL_TARGET_WRESTLER";
        case CAL_TARGET_CAYMAN:
            return "CAL_TARGET_CAYMAN";
        case CAL_TARGET_BARTS:
            return "CAL_TARGET_BARTS";
        case CAL_TARGET_RESERVED0:
            return "CAL_TARGET_RESERVED0";
        case CAL_TARGET_RESERVED1:
            return "CAL_TARGET_RESERVED1";
        case CAL_TARGET_RESERVED2:
            return "CAL_TARGET_RESERVED2";
        default:
            warn("Trying to show unhandled CALenumTarget: %d\n", x);
            return "Unknown CALtargetEnum";
    }
}

