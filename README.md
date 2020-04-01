# tviblindi

Xquartz necessary on macos.... so it seems


Witness complexes require CGAL library version 4.*.*

MacOS Catalina:
/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/sys/_types/_uuid_t.h:

28: /* #ifndef _UUID_T */

29: #define _UUID_T

30: #include <sys/_types.h> /* __darwin_uuid_t */

31: typedef __darwin_uuid_t	uuid_t;

32: /* #endif */ /* _UUID_T */
