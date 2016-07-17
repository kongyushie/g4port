#pragma once
#include "type.h"
#include "debug.h"
#include "element.h"
//#include "VCrossSectionDataSet.h"

typedef struct vector_dREAL {
	dREAL v[100];
	int length;
	int size;
}vector_dREAL;
extern void vector_dREAL_push(vector_dREAL* self,dREAL v);
extern int vector_dREAL_size(vector_dREAL* self);
extern void vector_dREAL_resize(vector_dREAL* self, int size);

typedef struct vector_VCrossSectionDataSet {
	dREAL v[100];
	int length;
	int size;
}vector_VCrossSectionDataSet;
/*
typedef struct vector_element_p {
	int v[100];
	int length;
	int size;
}vector_element_p;

typedef struct vector_isotope_p {
	int v[100];
	int length;
	int size;
}vector_isotope_p;
*/
#define SIZE 20
#define DEF_VECP(TYPE) \
typedef TYPE* TYPE##_p;\
DEF_VEC(TYPE##_p)

#define DEF_VEC(TYPE) \
typedef struct vector_##TYPE {\
	int length;\
		TYPE v[SIZE];\
} vector_##TYPE

#define VADD(ARR, item) \
if((ARR)->length!=SIZE) {\
		(ARR)->v[(ARR)->length++] = (item);\
} else {\
	exit(-1);\
}
