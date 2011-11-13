
#ifndef _CRLIBM_MSVC_HACKS_H_
#define _CRLIBM_MSVC_HACKS_H_

#ifdef _MSC_VER
typedef signed __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef signed __int64  int64_t;
typedef unsigned __int64 uint64_t;
#endif /* _MSC_VER */

#endif /* _CRLIBM_MSVC_HACKS_H_ */

