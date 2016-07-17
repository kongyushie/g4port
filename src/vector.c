#include "vector.h"

void vector_dREAL_add(vector_dREAL* self,dREAL value){
	if((*self).length < (*self).size+1){
		int index = (*self).size;
		(*self).v[index] = value;
		(*self).size++;
	}
	else {
		BLURT;
		printf("not dynamic initial \n");
		exit(0);	
	}
	return;
	
}

int vector_dREAL_size(vector_dREAL* self){
	return (*self).size;
}

void vector_dREAL_resize(vector_dREAL* self, int size){
	if(size < 100){
		(*self).length = size;
	}
	else {
		BLURT;
		printf("not dynamic initial \n");
		exit(0);
	}
	return;
}
