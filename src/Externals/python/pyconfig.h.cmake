/* pyconfig.h.in.  Generated from configure.in by autoheader.  */


#ifndef Py_PYCONFIG_H
#define Py_PYCONFIG_H


/* Define for AIX if your compiler is a genuine IBM xlC/xlC_r and you want
   support for AIX C++ shared extension modules. */
#cmakedefine AIX_GENUINE_CPLUSPLUS 1

/* Define this if you have AtheOS threads. */
#cmakedefine ATHEOS_THREADS 1

/* Define if you have the Mach cthreads package */
#cmakedefine C_THREADS 1

/* Define if C doubles are 64-bit IEEE 754 binary format, stored in ARM
   mixed-endian order (byte order 45670123) */
#cmakedefine DOUBLE_IS_ARM_MIXED_ENDIAN_IEEE754

/* Define if C doubles are 64-bit IEEE 754 binary format, stored with the most
   significant byte first */
#cmakedefine DOUBLE_IS_BIG_ENDIAN_IEEE754

/* Define if C doubles are 64-bit IEEE 754 binary format, stored with the
   least significant byte first */
#cmakedefine DOUBLE_IS_LITTLE_ENDIAN_IEEE754

/* Define if --enable-ipv6 is specified */
#cmakedefine ENABLE_IPV6 1

/* Define if getpgrp() must be called as getpgrp(0). */
#cmakedefine GETPGRP_HAVE_ARG 1

/* Define if gettimeofday() does not have second (timezone) argument This is
   the case on Motorola V4 (R40V4.2) */
#cmakedefine GETTIMEOFDAY_NO_TZ 1

/* Define to 1 if you have the `acosh' function. */
#cmakedefine HAVE_ACOSH 1

/* struct addrinfo (netdb.h) */
#cmakedefine HAVE_ADDRINFO 1

/* Define to 1 if you have the `alarm' function. */
#cmakedefine HAVE_ALARM 1

/* Define this if your time.h defines altzone. */
#cmakedefine HAVE_ALTZONE 1

/* Define to 1 if you have the `asinh' function. */
#cmakedefine HAVE_ASINH 1

/* Define to 1 if you have the <asm/types.h> header file. */
#cmakedefine HAVE_ASM_TYPES_H 1

/* Define to 1 if you have the `atanh' function. */
#cmakedefine HAVE_ATANH 1

/* Define if GCC supports __attribute__((format(PyArg_ParseTuple, 2, 3))) */
#cmakedefine HAVE_ATTRIBUTE_FORMAT_PARSETUPLE

/* Define to 1 if you have the `bind_textdomain_codeset' function. */
#cmakedefine HAVE_BIND_TEXTDOMAIN_CODESET 1

/* Define to 1 if you have the <bluetooth/bluetooth.h> header file. */
#cmakedefine HAVE_BLUETOOTH_BLUETOOTH_H 1

/* Define to 1 if you have the <bluetooth.h> header file. */
#cmakedefine HAVE_BLUETOOTH_H 1

/* Define if mbstowcs(NULL, "text", 0) does not return the number of wide
   chars that would be converted. */
#cmakedefine HAVE_BROKEN_MBSTOWCS

/* Define if nice() returns success/failure instead of the new priority. */
#cmakedefine HAVE_BROKEN_NICE 1

/* Define if poll() sets errno on invalid file descriptors. */
#cmakedefine HAVE_BROKEN_POLL 1

/* Define if the Posix semaphores do not work on your system */
#cmakedefine HAVE_BROKEN_POSIX_SEMAPHORES

/* Define if pthread_sigmask() does not work on your system. */
#cmakedefine HAVE_BROKEN_PTHREAD_SIGMASK

/* define to 1 if your sem_getvalue is broken. */
#cmakedefine HAVE_BROKEN_SEM_GETVALUE

/* Define this if you have the type _Bool. */
#cmakedefine HAVE_C99_BOOL 1

/* Define to 1 if you have the `chflags' function. */
#cmakedefine HAVE_CHFLAGS

/* Define to 1 if you have the `chown' function. */
#cmakedefine HAVE_CHOWN 1

/* Define if you have the 'chroot' function. */
#cmakedefine HAVE_CHROOT 1

/* Define to 1 if you have the `clock' function. */
#cmakedefine HAVE_CLOCK 1

/* Define to 1 if you have the `confstr' function. */
#cmakedefine HAVE_CONFSTR 1

/* Define to 1 if you have the <conio.h> header file. */
#cmakedefine HAVE_CONIO_H 1

/* Define to 1 if you have the `copysign' function. */
#cmakedefine HAVE_COPYSIGN

/* Define to 1 if you have the `ctermid' function. */
#cmakedefine HAVE_CTERMID 1

/* Define if you have the 'ctermid_r' function. */
#cmakedefine HAVE_CTERMID_R 1

/* Define to 1 if you have the <curses.h> header file. */
#cmakedefine HAVE_CURSES_H 1

/* Define if you have the 'is_term_resized' function. */
#cmakedefine HAVE_CURSES_IS_TERM_RESIZED 1

/* Define if you have the 'resizeterm' function. */
#cmakedefine HAVE_CURSES_RESIZETERM 1

/* Define if you have the 'resize_term' function. */
#cmakedefine HAVE_CURSES_RESIZE_TERM 1

/* Define to 1 if you have the declaration of `isfinite', and to 0 if you
   don't. */
#cmakedefine HAVE_DECL_ISFINITE

/* Define to 1 if you have the declaration of `isinf', and to 0 if you don't.
   */
#cmakedefine HAVE_DECL_ISINF

/* Define to 1 if you have the declaration of `isnan', and to 0 if you don't.
   */
#cmakedefine HAVE_DECL_ISNAN

/* Define to 1 if you have the declaration of `tzname', and to 0 if you don't.
   */
#cmakedefine HAVE_DECL_TZNAME

/* Define to 1 if you have the device macros. */
#cmakedefine HAVE_DEVICE_MACROS

/* Define if we have /dev/ptc. */
#cmakedefine HAVE_DEV_PTC

/* Define if we have /dev/ptmx. */
#cmakedefine HAVE_DEV_PTMX

/* Define to 1 if you have the <direct.h> header file. */
#cmakedefine HAVE_DIRECT_H

/* Define to 1 if you have the <dirent.h> header file, and it defines `DIR'.
   */
#cmakedefine HAVE_DIRENT_H

/* Define to 1 if you have the <dlfcn.h> header file. */
#cmakedefine HAVE_DLFCN_H

/* Define to 1 if you have the `dlopen' function. */
#cmakedefine HAVE_DLOPEN

/* Define to 1 if you have the `dup2' function. */
#cmakedefine HAVE_DUP2

/* Defined when any dynamic module loading is enabled. */
#undef HAVE_DYNAMIC_LOADING

/* Define if you have the 'epoll' functions. */
#cmakedefine HAVE_EPOLL

/* Define to 1 if you have the <errno.h> header file. */
#cmakedefine HAVE_ERRNO_H

/* Define to 1 if you have the `execv' function. */
#cmakedefine HAVE_EXECV

/* Define to 1 if you have the `expm1' function. */
#cmakedefine HAVE_EXPM1

/* Define if you have the 'fchdir' function. */
#cmakedefine HAVE_FCHDIR

/* Define to 1 if you have the `fchmod' function. */
#cmakedefine HAVE_FCHMOD

/* Define to 1 if you have the `fchown' function. */
#cmakedefine HAVE_FCHOWN

/* Define to 1 if you have the <fcntl.h> header file. */
#cmakedefine HAVE_FCNTL_H

/* Define if you have the 'fdatasync' function. */
#cmakedefine HAVE_FDATASYNC

/* Define to 1 if you have the `finite' function. */
#cmakedefine HAVE_FINITE

/* Define if you have the 'flock' function. */
#cmakedefine HAVE_FLOCK

/* Define to 1 if you have the `fork' function. */
#cmakedefine HAVE_FORK

/* Define to 1 if you have the `forkpty' function. */
#cmakedefine HAVE_FORKPTY

/* Define to 1 if you have the `fpathconf' function. */
#cmakedefine HAVE_FPATHCONF

/* Define to 1 if you have the `fseek64' function. */
#cmakedefine HAVE_FSEEK64

/* Define to 1 if you have the `fseeko' function. */
#cmakedefine HAVE_FSEEKO

/* Define to 1 if you have the `fstatvfs' function. */
#cmakedefine HAVE_FSTATVFS

/* Define if you have the 'fsync' function. */
#cmakedefine HAVE_FSYNC

/* Define to 1 if you have the `ftell64' function. */
#cmakedefine HAVE_FTELL64

/* Define to 1 if you have the `ftello' function. */
#cmakedefine HAVE_FTELLO

/* Define to 1 if you have the `ftime' function. */
#cmakedefine HAVE_FTIME

/* Define to 1 if you have the `ftruncate' function. */
#cmakedefine HAVE_FTRUNCATE

/* Define to 1 if you have the `gai_strerror' function. */
#cmakedefine HAVE_GAI_STRERROR

/* Define if we can use gcc inline assembler to get and set x87 control word
   */
#cmakedefine HAVE_GCC_ASM_FOR_X87

/* Define if you have the getaddrinfo function. */
#cmakedefine HAVE_GETADDRINFO

/* Define to 1 if you have the `getcwd' function. */
#cmakedefine HAVE_GETCWD

/* Define this if you have flockfile(), getc_unlocked(), and funlockfile() */
#cmakedefine HAVE_GETC_UNLOCKED

/* Define to 1 if you have the `getgroups' function. */
#cmakedefine HAVE_GETGROUPS

/* Define to 1 if you have the `gethostbyname' function. */
#cmakedefine HAVE_GETHOSTBYNAME

/* Define this if you have some version of gethostbyname_r() */
#cmakedefine HAVE_GETHOSTBYNAME_R

/* Define this if you have the 3-arg version of gethostbyname_r(). */
#cmakedefine HAVE_GETHOSTBYNAME_R_3_ARG

/* Define this if you have the 5-arg version of gethostbyname_r(). */
#cmakedefine HAVE_GETHOSTBYNAME_R_5_ARG

/* Define this if you have the 6-arg version of gethostbyname_r(). */
#cmakedefine HAVE_GETHOSTBYNAME_R_6_ARG

/* Define to 1 if you have the `getitimer' function. */
#cmakedefine HAVE_GETITIMER

/* Define to 1 if you have the `getloadavg' function. */
#cmakedefine HAVE_GETLOADAVG

/* Define to 1 if you have the `getlogin' function. */
#cmakedefine HAVE_GETLOGIN

/* Define to 1 if you have the `getnameinfo' function. */
#cmakedefine HAVE_GETNAMEINFO

/* Define if you have the 'getpagesize' function. */
#cmakedefine HAVE_GETPAGESIZE

/* Define to 1 if you have the `getpeername' function. */
#cmakedefine HAVE_GETPEERNAME

/* Define to 1 if you have the `getpgid' function. */
#cmakedefine HAVE_GETPGID

/* Define to 1 if you have the `getpgrp' function. */
#cmakedefine HAVE_GETPGRP

/* Define to 1 if you have the `getpid' function. */
#cmakedefine HAVE_GETPID

/* Define to 1 if you have the `getpriority' function. */
#cmakedefine HAVE_GETPRIORITY

/* Define to 1 if you have the `getpwent' function. */
#cmakedefine HAVE_GETPWENT

/* Define to 1 if you have the `getsid' function. */
#cmakedefine HAVE_GETSID

/* Define to 1 if you have the `getspent' function. */
#cmakedefine HAVE_GETSPENT

/* Define to 1 if you have the `getspnam' function. */
#cmakedefine HAVE_GETSPNAM

/* Define to 1 if you have the `gettimeofday' function. */
#cmakedefine HAVE_GETTIMEOFDAY

/* Define to 1 if you have the `getwd' function. */
#cmakedefine HAVE_GETWD

/* Define to 1 if you have the <grp.h> header file. */
#cmakedefine HAVE_GRP_H

/* Define if you have the 'hstrerror' function. */
#cmakedefine HAVE_HSTRERROR

/* Define to 1 if you have the `hypot' function. */
#cmakedefine HAVE_HYPOT

/* Define to 1 if you have the <ieeefp.h> header file. */
#cmakedefine HAVE_IEEEFP_H

/* Define if you have the 'inet_aton' function. */
#cmakedefine HAVE_INET_ATON

/* Define if you have the 'inet_pton' function. */
#cmakedefine HAVE_INET_PTON

/* Define to 1 if you have the <inttypes.h> header file. */
#cmakedefine HAVE_INTTYPES_H

/* Define to 1 if you have the <io.h> header file. */
#cmakedefine HAVE_IO_H

/* Define to 1 if you have the `kill' function. */
#cmakedefine HAVE_KILL

/* Define to 1 if you have the `killpg' function. */
#cmakedefine HAVE_KILLPG

/* Define if you have the 'kqueue' functions. */
#cmakedefine HAVE_KQUEUE

/* Define to 1 if you have the <langinfo.h> header file. */
#cmakedefine HAVE_LANGINFO_H

/* Defined to enable large file support when an off_t is bigger than a long
   and long long is available and at least as big as an off_t. You may need to
   add some flags for configuration and compilation to enable this mode. (For
   Solaris and Linux, the necessary defines are already defined.) */
#cmakedefine HAVE_LARGEFILE_SUPPORT @HAVE_LARGEFILE_SUPPORT@

/* Define to 1 if you have the `lchflags' function. */
#cmakedefine HAVE_LCHFLAGS

/* Define to 1 if you have the `lchmod' function. */
#cmakedefine HAVE_LCHMOD

/* Define to 1 if you have the `lchown' function. */
#cmakedefine HAVE_LCHOWN

/* Define to 1 if you have the `dl' library (-ldl). */
#cmakedefine HAVE_LIBDL

/* Define to 1 if you have the `dld' library (-ldld). */
#cmakedefine HAVE_LIBDLD

/* Define to 1 if you have the `ieee' library (-lieee). */
#cmakedefine HAVE_LIBIEEE

/* Define to 1 if you have the <libintl.h> header file. */
#cmakedefine HAVE_LIBINTL_H

/* Define if you have the readline library (-lreadline). */
#cmakedefine HAVE_LIBREADLINE

/* Define to 1 if you have the `resolv' library (-lresolv). */
#cmakedefine HAVE_LIBRESOLV

/* Define to 1 if you have the <libutil.h> header file. */
#cmakedefine HAVE_LIBUTIL_H

/* Define if you have the 'link' function. */
#cmakedefine HAVE_LINK

/* Define to 1 if you have the <linux/netlink.h> header file. */
#cmakedefine HAVE_LINUX_NETLINK_H

/* Define to 1 if you have the <linux/tipc.h> header file. */
#cmakedefine HAVE_LINUX_TIPC_H

/* Define to 1 if you have the `log1p' function. */
#cmakedefine HAVE_LOG1P

/* Define this if you have the type long double. */
#cmakedefine HAVE_LONG_DOUBLE

/* Define this if you have the type long long. */
#cmakedefine HAVE_LONG_LONG

/* Define to 1 if you have the `lstat' function. */
#cmakedefine HAVE_LSTAT

/* Define this if you have the makedev macro. */
#cmakedefine HAVE_MAKEDEV

/* Define to 1 if you have the `mbrtowc' function. */
#cmakedefine HAVE_MBRTOWC

/* Define to 1 if you have the `memmove' function. */
#cmakedefine HAVE_MEMMOVE

/* Define to 1 if you have the <memory.h> header file. */
#cmakedefine HAVE_MEMORY_H

/* Define to 1 if you have the `mkfifo' function. */
#cmakedefine HAVE_MKFIFO

/* Define to 1 if you have the `mknod' function. */
#cmakedefine HAVE_MKNOD

/* Define to 1 if you have the `mktime' function. */
#cmakedefine HAVE_MKTIME

/* Define to 1 if you have the `mremap' function. */
#cmakedefine HAVE_MREMAP

/* Define to 1 if you have the <ncurses.h> header file. */
#cmakedefine HAVE_NCURSES_H

/* Define to 1 if you have the <ndir.h> header file, and it defines `DIR'. */
#cmakedefine HAVE_NDIR_H

/* Define to 1 if you have the <netpacket/packet.h> header file. */
#cmakedefine HAVE_NETPACKET_PACKET_H

/* Define to 1 if you have the `nice' function. */
#cmakedefine HAVE_NICE

/* Define to 1 if you have the `openpty' function. */
#cmakedefine HAVE_OPENPTY

/* Define if compiling using MacOS X 10.5 SDK or later. */
#cmakedefine HAVE_OSX105_SDK

/* Define to 1 if you have the `pathconf' function. */
#cmakedefine HAVE_PATHCONF

/* Define to 1 if you have the `pause' function. */
#cmakedefine HAVE_PAUSE

/* Define to 1 if you have the `plock' function. */
#cmakedefine HAVE_PLOCK

/* Define to 1 if you have the `poll' function. */
#cmakedefine HAVE_POLL

/* Define to 1 if you have the <poll.h> header file. */
#cmakedefine HAVE_POLL_H

/* Define to 1 if you have the <process.h> header file. */
#cmakedefine HAVE_PROCESS_H

/* Define if your compiler supports function prototype */
#cmakedefine HAVE_PROTOTYPES

/* Define if you have GNU PTH threads. */
#cmakedefine HAVE_PTH

/* Defined for Solaris 2.6 bug in pthread header. */
#cmakedefine HAVE_PTHREAD_DESTRUCTOR

/* Define to 1 if you have the <pthread.h> header file. */
#cmakedefine HAVE_PTHREAD_H

/* Define to 1 if you have the `pthread_init' function. */
#cmakedefine HAVE_PTHREAD_INIT

/* Define to 1 if you have the `pthread_sigmask' function. */
#cmakedefine HAVE_PTHREAD_SIGMASK

/* Define to 1 if you have the <pty.h> header file. */
#cmakedefine HAVE_PTY_H

/* Define to 1 if you have the `putenv' function. */
#cmakedefine HAVE_PUTENV

/* Define to 1 if you have the `readlink' function. */
#cmakedefine HAVE_READLINK

/* Define to 1 if you have the `realpath' function. */
#cmakedefine HAVE_REALPATH

/* Define if you have readline 2.1 */
#cmakedefine HAVE_RL_CALLBACK

/* Define if you can turn off readline's signal handling. */
#cmakedefine HAVE_RL_CATCH_SIGNAL

/* Define if you have readline 2.2 */
#cmakedefine HAVE_RL_COMPLETION_APPEND_CHARACTER

/* Define if you have readline 4.0 */
#cmakedefine HAVE_RL_COMPLETION_DISPLAY_MATCHES_HOOK

/* Define if you have readline 4.2 */
#cmakedefine HAVE_RL_COMPLETION_MATCHES

/* Define if you have rl_completion_suppress_append */
#cmakedefine HAVE_RL_COMPLETION_SUPPRESS_APPEND

/* Define if you have readline 4.0 */
#cmakedefine HAVE_RL_PRE_INPUT_HOOK

/* Define to 1 if you have the `round' function. */
#cmakedefine HAVE_ROUND

/* Define to 1 if you have the `select' function. */
#cmakedefine HAVE_SELECT

/* Define to 1 if you have the `sem_getvalue' function. */
#cmakedefine HAVE_SEM_GETVALUE

/* Define to 1 if you have the `sem_open' function. */
#cmakedefine HAVE_SEM_OPEN

/* Define to 1 if you have the `sem_timedwait' function. */
#cmakedefine HAVE_SEM_TIMEDWAIT

/* Define to 1 if you have the `sem_unlink' function. */
#cmakedefine HAVE_SEM_UNLINK

/* Define to 1 if you have the `setegid' function. */
#cmakedefine HAVE_SETEGID

/* Define to 1 if you have the `seteuid' function. */
#cmakedefine HAVE_SETEUID

/* Define to 1 if you have the `setgid' function. */
#cmakedefine HAVE_SETGID

/* Define if you have the 'setgroups' function. */
#cmakedefine HAVE_SETGROUPS

/* Define to 1 if you have the `setitimer' function. */
#cmakedefine HAVE_SETITIMER

/* Define to 1 if you have the `setlocale' function. */
#cmakedefine HAVE_SETLOCALE

/* Define to 1 if you have the `setpgid' function. */
#cmakedefine HAVE_SETPGID

/* Define to 1 if you have the `setpgrp' function. */
#cmakedefine HAVE_SETPGRP

/* Define to 1 if you have the `setregid' function. */
#cmakedefine HAVE_SETREGID

/* Define to 1 if you have the `setreuid' function. */
#cmakedefine HAVE_SETREUID

/* Define to 1 if you have the `setsid' function. */
#cmakedefine HAVE_SETSID

/* Define to 1 if you have the `setuid' function. */
#cmakedefine HAVE_SETUID

/* Define to 1 if you have the `setvbuf' function. */
#cmakedefine HAVE_SETVBUF

/* Define to 1 if you have the <shadow.h> header file. */
#cmakedefine HAVE_SHADOW_H

/* Define to 1 if you have the `sigaction' function. */
#cmakedefine HAVE_SIGACTION

/* Define to 1 if you have the `siginterrupt' function. */
#cmakedefine HAVE_SIGINTERRUPT

/* Define to 1 if you have the <signal.h> header file. */
#cmakedefine HAVE_SIGNAL_H

/* Define to 1 if you have the `sigrelse' function. */
#cmakedefine HAVE_SIGRELSE

/* Define to 1 if you have the `snprintf' function. */
#cmakedefine HAVE_SNPRINTF

/* Define if sockaddr has sa_len member */
#cmakedefine HAVE_SOCKADDR_SA_LEN

/* struct sockaddr_storage (sys/socket.h) */
#cmakedefine HAVE_SOCKADDR_STORAGE

/* Define if you have the 'socketpair' function. */
#cmakedefine HAVE_SOCKETPAIR

/* Define if your compiler provides ssize_t */
#cmakedefine HAVE_SSIZE_T

/* Define to 1 if you have the `statvfs' function. */
#cmakedefine HAVE_STATVFS

/* Define if you have struct stat.st_mtim.tv_nsec */
#cmakedefine HAVE_STAT_TV_NSEC

/* Define if you have struct stat.st_mtimensec */
#cmakedefine HAVE_STAT_TV_NSEC2

/* Define if your compiler supports variable length function prototypes (e.g.
   void fprintf(FILE *, char *, ...);) *and* <stdarg.h> */
#cmakedefine HAVE_STDARG_PROTOTYPES

/* Define to 1 if you have the <stdint.h> header file. */
#cmakedefine HAVE_STDINT_H

/* Define to 1 if you have the <stdlib.h> header file. */
#cmakedefine HAVE_STDLIB_H

/* Define to 1 if you have the `strdup' function. */
#cmakedefine HAVE_STRDUP

/* Define to 1 if you have the `strftime' function. */
#cmakedefine HAVE_STRFTIME

/* Define to 1 if you have the <strings.h> header file. */
#cmakedefine HAVE_STRINGS_H

/* Define to 1 if you have the <string.h> header file. */
#cmakedefine HAVE_STRING_H

/* Define to 1 if you have the `strlcpy' function. */
#cmakedefine HAVE_STRLCPY

/* Define to 1 if you have the <stropts.h> header file. */
#cmakedefine HAVE_STROPTS_H

/* Define to 1 if `st_birthtime' is member of `struct stat'. */
#cmakedefine HAVE_STRUCT_STAT_ST_BIRTHTIME

/* Define to 1 if `st_blksize' is member of `struct stat'. */
#cmakedefine HAVE_STRUCT_STAT_ST_BLKSIZE

/* Define to 1 if `st_blocks' is member of `struct stat'. */
#cmakedefine HAVE_STRUCT_STAT_ST_BLOCKS

/* Define to 1 if `st_flags' is member of `struct stat'. */
#cmakedefine HAVE_STRUCT_STAT_ST_FLAGS

/* Define to 1 if `st_gen' is member of `struct stat'. */
#cmakedefine HAVE_STRUCT_STAT_ST_GEN

/* Define to 1 if `st_rdev' is member of `struct stat'. */
#cmakedefine HAVE_STRUCT_STAT_ST_RDEV

/* Define to 1 if `tm_zone' is member of `struct tm'. */
#cmakedefine HAVE_STRUCT_TM_TM_ZONE

/* Define to 1 if your `struct stat' has `st_blocks'. Deprecated, use
   `HAVE_STRUCT_STAT_ST_BLOCKS' instead. */
#cmakedefine HAVE_ST_BLOCKS

/* Define if you have the 'symlink' function. */
#cmakedefine HAVE_SYMLINK

/* Define to 1 if you have the `sysconf' function. */
#cmakedefine HAVE_SYSCONF

/* Define to 1 if you have the <sysexits.h> header file. */
#cmakedefine HAVE_SYSEXITS_H

/* Define to 1 if you have the <sys/audioio.h> header file. */
#cmakedefine HAVE_SYS_AUDIOIO_H

/* Define to 1 if you have the <sys/bsdtty.h> header file. */
#cmakedefine HAVE_SYS_BSDTTY_H

/* Define to 1 if you have the <sys/dir.h> header file, and it defines `DIR'.
   */
#cmakedefine HAVE_SYS_DIR_H

/* Define to 1 if you have the <sys/epoll.h> header file. */
#cmakedefine HAVE_SYS_EPOLL_H

/* Define to 1 if you have the <sys/event.h> header file. */
#cmakedefine HAVE_SYS_EVENT_H

/* Define to 1 if you have the <sys/file.h> header file. */
#cmakedefine HAVE_SYS_FILE_H

/* Define to 1 if you have the <sys/loadavg.h> header file. */
#cmakedefine HAVE_SYS_LOADAVG_H

/* Define to 1 if you have the <sys/lock.h> header file. */
#cmakedefine HAVE_SYS_LOCK_H

/* Define to 1 if you have the <sys/mkdev.h> header file. */
#cmakedefine HAVE_SYS_MKDEV_H

/* Define to 1 if you have the <sys/modem.h> header file. */
#cmakedefine HAVE_SYS_MODEM_H

/* Define to 1 if you have the <sys/ndir.h> header file, and it defines `DIR'.
   */
#cmakedefine HAVE_SYS_NDIR_H

/* Define to 1 if you have the <sys/param.h> header file. */
#cmakedefine HAVE_SYS_PARAM_H

/* Define to 1 if you have the <sys/poll.h> header file. */
#cmakedefine HAVE_SYS_POLL_H

/* Define to 1 if you have the <sys/resource.h> header file. */
#cmakedefine HAVE_SYS_RESOURCE_H

/* Define to 1 if you have the <sys/select.h> header file. */
#cmakedefine HAVE_SYS_SELECT_H

/* Define to 1 if you have the <sys/socket.h> header file. */
#cmakedefine HAVE_SYS_SOCKET_H

/* Define to 1 if you have the <sys/statvfs.h> header file. */
#cmakedefine HAVE_SYS_STATVFS_H

/* Define to 1 if you have the <sys/stat.h> header file. */
#cmakedefine HAVE_SYS_STAT_H

/* Define to 1 if you have the <sys/termio.h> header file. */
#cmakedefine HAVE_SYS_TERMIO_H

/* Define to 1 if you have the <sys/times.h> header file. */
#cmakedefine HAVE_SYS_TIMES_H

/* Define to 1 if you have the <sys/time.h> header file. */
#cmakedefine HAVE_SYS_TIME_H

/* Define to 1 if you have the <sys/types.h> header file. */
#cmakedefine HAVE_SYS_TYPES_H

/* Define to 1 if you have the <sys/un.h> header file. */
#cmakedefine HAVE_SYS_UN_H

/* Define to 1 if you have the <sys/utsname.h> header file. */
#cmakedefine HAVE_SYS_UTSNAME_H

/* Define to 1 if you have the <sys/wait.h> header file. */
#cmakedefine HAVE_SYS_WAIT_H

/* Define to 1 if you have the `tcgetpgrp' function. */
#cmakedefine HAVE_TCGETPGRP

/* Define to 1 if you have the `tcsetpgrp' function. */
#cmakedefine HAVE_TCSETPGRP

/* Define to 1 if you have the `tempnam' function. */
#cmakedefine HAVE_TEMPNAM

/* Define to 1 if you have the <termios.h> header file. */
#cmakedefine HAVE_TERMIOS_H

/* Define to 1 if you have the <term.h> header file. */
#cmakedefine HAVE_TERM_H

/* Define to 1 if you have the <thread.h> header file. */
#cmakedefine HAVE_THREAD_H

/* Define to 1 if you have the `timegm' function. */
#cmakedefine HAVE_TIMEGM

/* Define to 1 if you have the `times' function. */
#cmakedefine HAVE_TIMES

/* Define to 1 if you have the `tmpfile' function. */
#cmakedefine HAVE_TMPFILE

/* Define to 1 if you have the `tmpnam' function. */
#cmakedefine HAVE_TMPNAM

/* Define to 1 if you have the `tmpnam_r' function. */
#cmakedefine HAVE_TMPNAM_R

/* Define to 1 if your `struct tm' has `tm_zone'. Deprecated, use
   `HAVE_STRUCT_TM_TM_ZONE' instead. */
#cmakedefine HAVE_TM_ZONE

/* Define to 1 if you have the `truncate' function. */
#cmakedefine HAVE_TRUNCATE

/* Define to 1 if you don't have `tm_zone' but do have the external array
   `tzname'. */
#cmakedefine HAVE_TZNAME

/* Define this if you have tcl and TCL_UTF_MAX==6 */
#cmakedefine HAVE_UCS4_TCL

/* Define to 1 if the system has the type `uintptr_t'. */
#cmakedefine HAVE_UINTPTR_T

/* Define to 1 if you have the `uname' function. */
#cmakedefine HAVE_UNAME

/* Define to 1 if you have the <unistd.h> header file. */
#cmakedefine HAVE_UNISTD_H

/* Define to 1 if you have the `unsetenv' function. */
#cmakedefine HAVE_UNSETENV

/* Define if you have a useable wchar_t type defined in wchar.h; useable means
   wchar_t must be an unsigned type with at least 16 bits. (see
   Include/unicodeobject.h). */
#cmakedefine HAVE_USABLE_WCHAR_T

/* Define to 1 if you have the `utimes' function. */
#cmakedefine HAVE_UTIMES

/* Define to 1 if you have the <utime.h> header file. */
#cmakedefine HAVE_UTIME_H

/* Define to 1 if you have the `wait3' function. */
#cmakedefine HAVE_WAIT3

/* Define to 1 if you have the `wait4' function. */
#cmakedefine HAVE_WAIT4

/* Define to 1 if you have the `waitpid' function. */
#cmakedefine HAVE_WAITPID

/* Define if the compiler provides a wchar.h header file. */
#cmakedefine HAVE_WCHAR_H

/* Define to 1 if you have the `wcscoll' function. */
#cmakedefine HAVE_WCSCOLL

/* Define to 1 if you have the `wcsftime' function. */
#cmakedefine HAVE_WCSFTIME

/* Define to 1 if you have the `wcsxfrm' function. */
#cmakedefine HAVE_WCSXFRM

/* Define if tzset() actually switches the local timezone in a meaningful way.
   */
#cmakedefine HAVE_WORKING_TZSET

/* Define if the zlib library has inflateCopy */
#cmakedefine HAVE_ZLIB_COPY

/* Define to 1 if you have the `_getpty' function. */
#cmakedefine HAVE__GETPTY

/* Define if you are using Mach cthreads directly under /include */
#cmakedefine HURD_C_THREADS

/* Define if you are using Mach cthreads under mach / */
#cmakedefine MACH_C_THREADS

/* Define to 1 if `major', `minor', and `makedev' are declared in <mkdev.h>.
   */
#cmakedefine MAJOR_IN_MKDEV

/* Define to 1 if `major', `minor', and `makedev' are declared in
   <sysmacros.h>. */
#cmakedefine MAJOR_IN_SYSMACROS

/* Define if mvwdelch in curses.h is an expression. */
#cmakedefine MVWDELCH_IS_EXPRESSION

/* Define to the address where bug reports for this package should be sent. */
#cmakedefine PACKAGE_BUGREPORT

/* Define to the full name of this package. */
#cmakedefine PACKAGE_NAME @PACKAGE_NAME@

/* Define to the full name and version of this package. */
#cmakedefine PACKAGE_STRING @PACKAGE_STRING@

/* Define to the one symbol short name of this package. */
#cmakedefine PACKAGE_TARNAME @PACKAGE_TARNAME@

/* Define to the version of this package. */
#cmakedefine PACKAGE_VERSION @PACKAGE_VERSION@

/* Define if POSIX semaphores aren't enabled on your system */
#cmakedefine POSIX_SEMAPHORES_NOT_ENABLED

/* Defined if PTHREAD_SCOPE_SYSTEM supported. */
#cmakedefine PTHREAD_SYSTEM_SCHED_SUPPORTED

/* Define as the preferred size in bits of long digits */
#cmakedefine PYLONG_BITS_IN_DIGIT

/* Define to printf format modifier for Py_ssize_t */
#cmakedefine PY_FORMAT_SIZE_T

/* Define as the integral type used for Unicode representation. */
#cmakedefine PY_UNICODE_TYPE @PY_UNICODE_TYPE@

/* Define if you want to build an interpreter with many run-time checks. */
#cmakedefine Py_DEBUG

/* Defined if Python is built as a shared library. */
#cmakedefine Py_ENABLE_SHARED

/* Define as the size of the unicode type. */
#cmakedefine Py_UNICODE_SIZE @Py_UNICODE_SIZE@

/* Define as the return type of signal handlers (`int' or `void'). */
#cmakedefine RETSIGTYPE @RETSIGTYPE@

/* Define if setpgrp() must be called as setpgrp(0, 0). */
#cmakedefine SETPGRP_HAVE_ARG

/* Define this to be extension of shared libraries (including the dot!). */
#cmakedefine SHLIB_EXT "@CMAKE_SHARED_LIBRARY_SUFFIX@"

/* Define if i>>j for signed int i does not extend the sign bit when i < 0 */
#cmakedefine SIGNED_RIGHT_SHIFT_ZERO_FILLS

/* The size of `double', as computed by sizeof. */
#cmakedefine SIZEOF_DOUBLE @SIZEOF_DOUBLE@

/* The size of `float', as computed by sizeof. */
#cmakedefine SIZEOF_FLOAT @SIZEOF_FLOAT@

/* The size of `fpos_t', as computed by sizeof. */
#cmakedefine SIZEOF_FPOS_T @SIZEOF_FPOS_T@

/* The size of `int', as computed by sizeof. */
#cmakedefine SIZEOF_INT @SIZEOF_INT@

/* The size of `long', as computed by sizeof. */
#cmakedefine SIZEOF_LONG @SIZEOF_LONG@

/* The size of `long double', as computed by sizeof. */
#cmakedefine SIZEOF_LONG_DOUBLE @SIZEOF_LONG_DOUBLE@

/* The size of `long long', as computed by sizeof. */
#cmakedefine SIZEOF_LONG_LONG @SIZEOF_LONG_LONG@

/* The number of bytes in an off_t. */
#cmakedefine SIZEOF_OFF_T @SIZEOF_OFF_T@

/* The size of `pid_t', as computed by sizeof. */
#cmakedefine SIZEOF_PID_T @SIZEOF_PID_T@

/* The number of bytes in a pthread_t. */
#cmakedefine SIZEOF_PTHREAD_T @SIZEOF_PTHREAD_T@

/* The size of `short', as computed by sizeof. */
#cmakedefine SIZEOF_SHORT @SIZEOF_SHORT@

/* The size of `size_t', as computed by sizeof. */
#cmakedefine SIZEOF_SIZE_T @SIZEOF_SIZE_T@

/* The number of bytes in a time_t. */
#cmakedefine SIZEOF_TIME_T @SIZEOF_TIME_T@

/* The size of `uintptr_t', as computed by sizeof. */
#cmakedefine SIZEOF_UINTPTR_T @SIZEOF_UINTPTR_T@

/* The size of `void *', as computed by sizeof. */
#cmakedefine SIZEOF_VOID_P @SIZEOF_VOID_P@

/* The size of `wchar_t', as computed by sizeof. */
#cmakedefine SIZEOF_WCHAR_T @SIZEOF_WCHAR_T@

/* The size of `_Bool', as computed by sizeof. */
#cmakedefine SIZEOF__BOOL @SIZEOF__BOOL@

/* Define to 1 if you have the ANSI C header files. */
#cmakedefine STDC_HEADERS

/* Define if you can safely include both <sys/select.h> and <sys/time.h>
   (which you can't on SCO ODT 3.0). */
#cmakedefine SYS_SELECT_WITH_SYS_TIME

/* Define if tanh(-0.) is -0., or if platform doesn't have signed zeros */
#cmakedefine TANH_PRESERVES_ZERO_SIGN

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#cmakedefine TIME_WITH_SYS_TIME

/* Define to 1 if your <sys/time.h> declares `struct tm'. */
#cmakedefine TM_IN_SYS_TIME

/* Define if you want to use computed gotos in ceval.c. */
#cmakedefine USE_COMPUTED_GOTOS

/* Define if a va_list is an array of some kind */
#cmakedefine VA_LIST_IS_ARRAY

/* Define if you want SIGFPE handled (see Include/pyfpe.h). */
#cmakedefine WANT_SIGFPE_HANDLER

/* Define if you want wctype.h functions to be used instead of the one
   supplied by Python itself. (see Include/unicodectype.h). */
#cmakedefine WANT_WCTYPE_FUNCTIONS

/* Define if WINDOW in curses.h offers a field _flags. */
#cmakedefine WINDOW_HAS_FLAGS

/* Define if you want documentation strings in extension modules */
#cmakedefine WITH_DOC_STRINGS

/* Define if you want to use the new-style (Openstep, Rhapsody, MacOS) dynamic
   linker (dyld) instead of the old-style (NextStep) dynamic linker (rld).
   Dyld is necessary to support frameworks. */
#cmakedefine WITH_DYLD

/* Define to 1 if libintl is needed for locale functions. */
#cmakedefine WITH_LIBINTL

/* Define if you want to produce an OpenStep/Rhapsody framework (shared
   library plus accessory files). */
#cmakedefine WITH_NEXT_FRAMEWORK

/* Define if you want to compile in Python-specific mallocs */
#cmakedefine WITH_PYMALLOC

/* Define if you want to compile in rudimentary thread support */
#cmakedefine WITH_THREAD

/* Define to profile with the Pentium timestamp counter */
#cmakedefine WITH_TSC

/* Define to 1 if your processor stores words with the most significant byte
   first (like Motorola and SPARC, unlike Intel and VAX). */
#cmakedefine WORDS_BIGENDIAN

/* Define if arithmetic is subject to x87-style double rounding issue */
#cmakedefine X87_DOUBLE_ROUNDING

/* Define to 1 if on AIX 3.
   System headers sometimes define this.
   We just want to avoid a redefinition error message.  */
#ifndef _ALL_SOURCE
# undef _ALL_SOURCE
#endif

/* Define on OpenBSD to activate all library features */
#cmakedefine _BSD_SOURCE

/* Define on Irix to enable u_int */
#cmakedefine _BSD_TYPES

/* Define on Darwin to activate all library features */
#cmakedefine _DARWIN_C_SOURCE

/* This must be set to 64 on some systems to enable large file support. */
#cmakedefine _FILE_OFFSET_BITS @_FILE_OFFSET_BITS@

/* Define on Linux to activate all library features */
#cmakedefine _GNU_SOURCE

/* This must be defined on some systems to enable large file support. */
#cmakedefine _LARGEFILE_SOURCE @_LARGEFILE_SOURCE@

/* Define on NetBSD to activate all library features */
#cmakedefine _NETBSD_SOURCE

/* Define _OSF_SOURCE to get the makedev macro. */
#cmakedefine _OSF_SOURCE

/* Define to activate features from IEEE Stds 1003.1-2001 */
#cmakedefine _POSIX_C_SOURCE @_POSIX_C_SOURCE@

/* Define if you have POSIX threads, and your system does not define that. */
#cmakedefine _POSIX_THREADS

/* Define to force use of thread-safe errno, h_errno, and other functions */
#cmakedefine _REENTRANT

/* Define for Solaris 2.5.1 so the uint32_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef was allowed, the
   #define below would cause a syntax error. */
#cmakedefine _UINT32_T

/* Define for Solaris 2.5.1 so the uint64_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef was allowed, the
   #define below would cause a syntax error. */
#cmakedefine _UINT64_T

/* Define to the level of X/Open that your system supports */
#cmakedefine _XOPEN_SOURCE @_XOPEN_SOURCE@

/* Define to activate Unix95-and-earlier features */
#cmakedefine _XOPEN_SOURCE_EXTENDED

/* Define on FreeBSD to activate all library features */
#cmakedefine __BSD_VISIBLE

/* Define to 1 if type `char' is unsigned and you are not using gcc.  */
#ifndef __CHAR_UNSIGNED__
# undef __CHAR_UNSIGNED__
#endif

/* Defined on Solaris to see additional function prototypes. */
#cmakedefine __EXTENSIONS__

/* Define to 'long' if <time.h> doesn't define. */
#cmakedefine clock_t @clock_t@

/* Define to empty if `const' does not conform to ANSI C. */
#cmakedefine const

/* Define to `int' if <sys/types.h> doesn't define. */
#cmakedefine gid_t @gid_t@

/* Define to the type of a signed integer type of width exactly 32 bits if
   such a type exists and the standard includes do not define it. */
#cmakedefine int32_t

/* Define to the type of a signed integer type of width exactly 64 bits if
   such a type exists and the standard includes do not define it. */
#cmakedefine int64_t

/* Define to `int' if <sys/types.h> does not define. */
#cmakedefine mode_t @mode_t@

/* Define to `long int' if <sys/types.h> does not define. */
#cmakedefine off_t @off_t@

/* Define to `int' if <sys/types.h> does not define. */
#cmakedefine pid_t @pid_t@

/* Define to empty if the keyword does not work. */
#cmakedefine signed

/* Define to `unsigned int' if <sys/types.h> does not define. */
#cmakedefine size_t @size_t@

/* Define to `int' if <sys/socket.h> does not define. */
#cmakedefine socklen_t @socklen_t@

/* Define to `int' if <sys/types.h> doesn't define. */
#cmakedefine uid_t @uid_t@

/* Define to the type of an unsigned integer type of width exactly 32 bits if
   such a type exists and the standard includes do not define it. */
#cmakedefine uint32_t

/* Define to the type of an unsigned integer type of width exactly 64 bits if
   such a type exists and the standard includes do not define it. */
#cmakedefine uint64_t

/* Define to empty if the keyword does not work. */
#cmakedefine volatile


/* Define the macros needed if on a UnixWare 7.x system. */
#if defined(__USLC__) && defined(__SCO_VERSION__)
#define STRICT_SYSV_CURSES /* Don't use ncurses extensions */
#endif

#endif /*Py_PYCONFIG_H*/

