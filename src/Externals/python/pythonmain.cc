/* Minimal main program -- everything is loaded from the library */

#include "Python.h"
#include <locale.h>
#include <string>

#ifdef __FreeBSD__
#include <floatingpoint.h>
#endif

#ifdef MS_WINDOWS
int
wmain(int argc, wchar_t **argv)
{
	std::wstring module_path( argv[ 0 ] );
	std::wstring::size_type pos = module_path.find_last_of( L'/' );
	module_path = module_path.substr( 0, pos + 1 );
	module_path += PYTHONPATH;
	Py_SetPath( module_path.c_str() );
	return Py_Main(argc, argv);
}
#else

#ifdef __APPLE__
extern "C"
{
	extern wchar_t* _Py_DecodeUTF8_surrogateescape(const char *s, Py_ssize_t size);
}
#endif

int
main(int argc, char **argv)
{
    wchar_t **argv_copy = (wchar_t **)PyMem_Malloc(sizeof(wchar_t*)*argc);
    /* We need a second copies, as Python might modify the first one. */
    wchar_t **argv_copy2 = (wchar_t **)PyMem_Malloc(sizeof(wchar_t*)*argc);
	int i, res;
    char *oldloc;
    /* 754 requires that FP exceptions run in "no stop" mode by default,
     * and until C vendors implement C99's ways to control FP exceptions,
     * Python requires non-stop mode.  Alas, some platforms enable FP
     * exceptions by default.  Here we disable them.
     */
#ifdef __FreeBSD__
    fp_except_t m;

    m = fpgetmask();
    fpsetmask(m & ~FP_X_OFL);
#endif
    if (!argv_copy || !argv_copy2) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }
    oldloc = strdup(setlocale(LC_ALL, NULL));
    setlocale(LC_ALL, "");
    for (i = 0; i < argc; i++) {
#ifdef __APPLE__
        argv_copy[i] = _Py_DecodeUTF8_surrogateescape(argv[i], strlen(argv[i]));
#else
        argv_copy[i] = _Py_char2wchar(argv[i], NULL);
#endif
        if (!argv_copy[i])
            return 1;
        argv_copy2[i] = argv_copy[i];
    }
    setlocale(LC_ALL, oldloc);
    free(oldloc);
	
	// Figure out the module search path
	std::wstring module_path( argv_copy[ 0 ] );
	std::wstring::size_type pos = module_path.find_last_of( L'/' );
	module_path = module_path.substr( 0, pos + 1 );
	module_path += PYTHONPATH;
	Py_IgnoreEnvironmentFlag = 1;
	Py_NoSiteFlag = 1;
	Py_InspectFlag = 1;
	Py_SetPath( module_path.c_str() );
    res = Py_Main(argc, argv_copy);
    for (i = 0; i < argc; i++) {
        PyMem_Free(argv_copy2[i]);
    }
    PyMem_Free(argv_copy);
    PyMem_Free(argv_copy2);
    return res;
}
#endif
