module m_log

    contains

#ifdef UNIX
include "./stdout_unix.inc"
#else
include "./stdout_win32.inc"
#endif

end module m_log