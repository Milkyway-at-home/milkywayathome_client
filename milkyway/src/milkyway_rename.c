/*
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "milkyway_config.h"

#if HAVE_SYS_TIME_H
  #include <sys/time.h>
#endif

#if HAVE_SYS_STAT_H
  #include <sys/stat.h>
#endif

#include "milkyway_util.h"
#include "milkyway_rename.h"
#include "milkyway_boinc_util.h"




#ifdef _WIN32

static int transactionFuncsInit = FALSE;
static int transactionFuncsOK = FALSE;

static HANDLE (WINAPI *__CreateTransaction) (LPSECURITY_ATTRIBUTES,
                                             LPGUID,
                                             DWORD,
                                             DWORD,
                                             DWORD,
                                             DWORD,
                                             LPWSTR) = NULL;

static BOOL (WINAPI *__MoveFileTransacted) (LPCTSTR,
                                            LPCTSTR,
                                            LPPROGRESS_ROUTINE,
                                            LPVOID,
                                            DWORD,
                                            HANDLE) = NULL;

static BOOL (WINAPI *__CommitTransaction) (HANDLE TransactionHandle) = NULL;



/* The transactional stuff is only available on Vista and later */
static void initW32TransactionalFunctions()
{
    HMODULE ktm32Lib;
    HMODULE kernel32Lib;

    transactionFuncsInit = TRUE;

    kernel32Lib = LoadLibrary("Kernel32.dll");
    if (!kernel32Lib)
    {
        mwPerrorW32("Could not load Kernel32.dll");
        return;
    }

    ktm32Lib = LoadLibrary("KtmW32.dll");
    if (!ktm32Lib)
    {
        mwPerrorW32("Could not load Ktm32.dll");
        return;
    }

    __CreateTransaction = GetProcAddress(ktm32Lib, "CreateTransaction");
    __CommitTransaction = GetProcAddress(ktm32Lib, "CommitTransaction");
    __MoveFileTransacted = GetProcAddress(kernel32Lib, "MoveFileTransactedA");

    transactionFuncsOK = (__CreateTransaction && __MoveFileTransacted && __CommitTransaction);

    if (!transactionFuncsOK)
    {
        mw_printf("Failed to get transaction functions\n");
    }
}

static int mw_rename_w32_fallback(const char* oldf, const char* newf)
{
    if (MoveFileExA(oldf, newf, MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
        return 0;
    return GetLastError();
}

static int mw_rename_w32_atomic(const char* oldf, const char* newf)
{
    HANDLE tx;

    tx = __CreateTransaction(NULL, NULL, 0, 0, 0, 0, L"AtomicFileRenameTransaction");
    if (!tx)
    {
        mwPerrorW32("Failed to create transaction for renaming '%s' to '%s'", oldf, newf);
        return 1;
    }

    if (!__MoveFileTransacted(oldf, newf, NULL, NULL, MOVEFILE_REPLACE_EXISTING, tx))
    {
        mwPerrorW32("Failed to move file '%s' to '%s'", oldf, newf);
        return 1;
    }

    if (!__CommitTransaction(tx))
    {
        mwPerrorW32("Failed to commit move of '%s' to '%s'", oldf, newf);
        return 1;
    }

    if (!CloseHandle(tx))
    {
        mwPerrorW32("Failed to close transaction handle for move of '%s' to '%s'", oldf, newf);
        return 1;
    }

    return 0;
}

static int mw_rename_w32(const char* oldf, const char* newf)
{
    /* It turns out that rename() does exist although it doesn't behave
    properly and errors if the destination file already exists which is
    wrong. This isn't quite atomic like it's supposed to be. */

    if (!transactionFuncsInit)
    {
        initW32TransactionalFunctions();
    }

    if (transactionFuncsOK)
    {
        return mw_rename_w32_atomic(oldf, newf);
    }
    else
    {
        return mw_rename_w32_fallback(oldf, newf);
    }
}
#endif /* _WIN32 */


/* Temporary stuff until patch BOINC */
#if BOINC_APPLICATION
static int mw_boinc_rename_aux(const char* oldf, const char* newf)
{
#ifdef _WIN32
    return mw_rename_w32(oldf, newf);
#else
    return rename(oldf, newf);
#endif
}

/* Pretty much boinc_rename() mangled a bit to fit here temporarily */
static int mw_boinc_rename(const char* old, const char* newf)
{
    int retval;
    const double fileRetryInterval = 5;

    retval = mw_boinc_rename_aux(old, newf);
    if (retval)
    {
        double start = mwGetTime();
        do
        {
            mw_boinc_sleep(2.0 * (double) rand() / (double) RAND_MAX);
            retval = mw_boinc_rename_aux(old, newf);
            if (!retval)
                break;
        }
        while (mwGetTime() < start + fileRetryInterval);
    }

    return retval;
}

int mw_rename(const char* oldf, const char* newf)
{
    return mw_boinc_rename(oldf, newf);
}

#else /* !BOINC_APPLICATION */

int mw_rename(const char* oldf, const char* newf)
{
  #ifndef _WIN32
    return rename(oldf, newf);
  #else
    return mw_rename_w32(oldf, newf);
  #endif /* _WIN32 */
}

#endif /* BOINC_APPLICATON */


