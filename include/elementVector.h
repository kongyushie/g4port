#pragma once
#include "element.h"
typedef struct vector_element_p {
	element* v[100];
	int length;
	int size;
}vector_element_p;
//vector<element*>
typedef vector_element_p elementVector;



