#ifndef _DEBUG_H_
#define _DEBUG_H_

/*
 * debug.h
 * macros for debugging.
 *
 * Created 30 Jan. 2007, Copyright (C) Robert R. Snapp
 */

#if (DEBUG > 0)
//#define debug(exp) exp
#define debug(n,exp) {if (n <= DEBUG) { exp; }}
#else
//#define debug(exp)
#define debug(n,exp)
#endif



#endif /* _DEBUG_H_ */
