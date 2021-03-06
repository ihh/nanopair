#ifndef STRINGMAP_INCLUDED
#define STRINGMAP_INCLUDED

#include "rbtree.h"
#include "list.h"
#include "vector.h"

/* Strings and string containers.
   A "string" object is just a C string of type char*.
*/

/* container functions for strings */
void* StringNew(const char *a);
void* StringCopy(void* a);
void StringDelete(void* a);
int StringCompare(void* a, void* b);
void StringPrint(FILE*, void* a);

/* more string utilities */
char* StringConcat (const char *a, const char *b);

/* mappings from strings to values - basically wrappers for RBTree functions */
typedef RBTree StringMap;
typedef RBNode StringMapNode;
#define newStringMap(ValueCopyFunc,ValueDestroyFunc,ValuePrintFunc) ((StringMap*) newRBTree (StringCompare, StringCopy, ValueCopyFunc, StringDelete, ValueDestroyFunc, StringPrint, ValuePrintFunc))
#define deleteStringMap(STRINGMAPPTR) deleteRBTree ((RBTree*) STRINGMAPPTR)
#define StringMapSet(STRINGMAPPTR,STRING,VALUE) ((StringMapNode*) RBTreeSet ((RBTree*) STRINGMAPPTR, (void*) StringNew(STRING), VALUE))
#define StringMapInsert(STRINGMAPPTR,STRING,VALUE) ((StringMapNode*) RBTreeInsert ((RBTree*) STRINGMAPPTR, (void*) StringNew(STRING), VALUE))
#define StringMapErase(STRINGMAPPTR,STRING) RBTreeErase ((RBTree*) STRINGMAPPTR, (void*) STRING)
#define StringMapFind(STRINGMAPPTR,STRING) ((StringMapNode*) RBTreeFind ((RBTree*) STRINGMAPPTR, (void*) STRING))

/* typedefs & macros for StringSet, a value-less StringMap */
typedef StringMap StringSet;
typedef StringMapNode StringSetNode;
#define newStringSet() ((StringSet*) newStringMap (NullCopyFunction, NullDestroyFunction, NullPrintFunction))
#define deleteStringSet(STRINGSETPTR) deleteStringMap((StringMap*)STRINGSETPTR)
#define StringSetInsert(STRINGSETPTR,STRING) ((StringSetNode*) StringMapInsert((StringMap*)STRINGSETPTR,STRING,NULL))
#define StringSetErase(STRINGSETPTR,STRING) StringMapErase((StringMap*)STRINGSETPTR,STRING)
#define StringSetFind(STRINGSETPTR,STRING) ((StringSetNode*) StringMapFind((StringMap*)STRINGSETPTR,STRING))

/* typedefs & macros for Dictionary, a map from Strings to Strings */
typedef StringMap Dictionary;
typedef StringMapNode DictionaryNode;
#define newDictionary() ((Dictionary*) newStringMap (StringCopy, StringDelete, StringPrint))
#define deleteDictionary(DICTPTR) deleteStringMap((StringMap*)DICTPTR)
#define DictionarySet(DICTPTR,STRING1,STRING2) ((DictionaryNode*) StringMapSet((StringMap*)DICTPTR,STRING1,STRING2))
#define DictionaryInsert(DICTPTR,STRING1,STRING2) ((DictionaryNode*) StringMapInsert((StringMap*)DICTPTR,STRING1,STRING2))
#define DictionaryErase(DICTPTR,STRING) StringMapErase((StringMap*)DICTPTR,STRING)
#define DictionaryFind(DICTPTR,STRING) ((DictionaryNode*) StringMapFind((StringMap*)DICTPTR,STRING))

/* typedefs & macros for StringIntMap, a map from Strings to integers */
typedef StringMap StringIntMap;
typedef StringMapNode StringIntMapNode;
#define newStringIntMap() ((StringIntMap*) newStringMap (IntCopy, IntDelete, IntPrint))
#define deleteStringIntMap(SIMPTR) deleteStringMap((StringMap*)SIMPTR)
#define StringIntMapSet(SIMPTR,STRING,INT) ((StringIntMapNode*) StringMapSet((StringMap*)SIMPTR,STRING,IntNew(INT)))
#define StringIntMapInsert(SIMPTR,STRING,INT) ((StringIntMapNode*) StringMapInsert((StringMap*)SIMPTR,STRING,IntNew(INT)))
#define StringIntMapErase(SIMPTR,STRING) StringMapErase((StringMap*)SIMPTR,STRING)
#define StringIntMapFind(SIMPTR,STRING) ((StringIntMapNode*) StringMapFind((StringMap*)SIMPTR,STRING))

/* typedefs & macros for StringDoubleMap, a map from Strings to doubles */
typedef StringMap StringDoubleMap;
typedef StringMapNode StringDoubleMapNode;
#define newStringDoubleMap() ((StringDoubleMap*) newStringMap (DoubleCopy, DoubleDelete, DoublePrint))
#define deleteStringDoubleMap(SIMPTR) deleteStringMap((StringMap*)SIMPTR)
#define StringDoubleMapSet(SIMPTR,STRING,DBL) ((StringDoubleMapNode*) StringMapSet((StringMap*)SIMPTR,STRING,DoubleNew(DBL)))
#define StringDoubleMapInsert(SIMPTR,STRING,DBL) ((StringDoubleMapNode*) StringMapInsert((StringMap*)SIMPTR,STRING,DoubleNew(DBL)))
#define StringDoubleMapErase(SIMPTR,STRING) StringMapErase((StringMap*)SIMPTR,STRING)
#define StringDoubleMapFind(SIMPTR,STRING) ((StringDoubleMapNode*) StringMapFind((StringMap*)SIMPTR,STRING))

/* StringVector */
typedef Vector StringVector;
#define newStringVector() ((StringVector*) newVector (StringCopy, StringDelete, StringPrint))
#define deleteStringVector(STRINGVECPTR) deleteVector ((Vector*) STRINGVECPTR)
#define StringVectorGet(STRINGVECPTR,N) ((const char*) VectorGet ((Vector*) STRINGVECPTR, N))
#define StringVectorSet(STRINGVECPTR,N,STRING) VectorSet ((Vector*) STRINGVECPTR, N, (void*) StringNew(STRING))
#define StringVectorReserve(STRINGVECPTR,N) VectorReserve ((Vector*) STRINGVECPTR, N)
#define StringVectorPushBack(STRINGVECPTR,STRING) VectorPushBack ((Vector*) STRINGVECPTR, (void*) StringNew(STRING))
#define StringVectorSize(STRINGVECPTR) VectorSize ((Vector*) STRINGVECPTR)


#endif /* STRINGMAP_INCLUDED */
