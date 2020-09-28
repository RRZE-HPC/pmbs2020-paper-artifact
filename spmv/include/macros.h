#ifndef _MACROS_H
#define _MACROS_H

/*
 * Concatenate preprocessor tokens A and B without expanding macro definitions
 * (however, if invoked from a macro, macro arguments are expanded).
 */
#define PPCAT_NX(A, B) A ## B

/*
 * Concatenate preprocessor tokens A and B after macro-expanding them.
 */
#define PPCAT(A, B) PPCAT_NX(A, B)

/*
 * Turn A into a string literal without expanding macro definitions
 * (however, if invoked from a macro, macro arguments are expanded).
 */
#define STRINGIZE_NX(A) #A

/*
 * Turn A into a string literal after macro-expanding it.
 */
#define STRINGIZE(A) STRINGIZE_NX(A)

#ifdef USE_SCC

#define BEGIN_SCC(_arr_)\
	_Pragma("statement scache_isolate_way L2=10 L1=2")\
	_Pragma(STRINGIZE(statement scache_isolate_assign _arr_))\
	;\

#define END_SCC\
	_Pragma("statement end_scache_isolate_assign")\
	_Pragma("statement end_scache_isolate_way")\


#else

#define BEGIN_SCC(_arr_)


#define END_SCC

#endif

#endif
