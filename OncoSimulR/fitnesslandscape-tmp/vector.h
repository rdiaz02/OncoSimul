#ifndef VECTOR_H
#define VECTOR_H

typedef struct vector {
    int **data;
	int *length;
    size_t capacity;
    size_t size;
} vector;

void vector_init(vector *v, size_t initialCapacity);
int vector_total(vector *v);
void vector_resize(vector *v, size_t capacity);
int vector_lengthElement(vector *v, int index);
void vector_add(vector * v, int * data, int length);
void vector_set(vector * v, int, int * data, int length);
int *vector_get(vector * v, int);
void vector_delete(vector * v, int);
void vector_free(vector * v);

#endif
