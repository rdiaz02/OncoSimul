#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vector.h"

void vector_init(vector *v, size_t initialCapacity)
{
    v->capacity = initialCapacity;
    v->size = 0;
	v->data = (int **)malloc(initialCapacity * sizeof(int*));
	v->length = (int *)malloc(initialCapacity * sizeof(int));
}

int vector_total(vector *v)
{
    return v->size;
}

int vector_lengthElement(vector *v, int index)
{
	return v->length[index];
}


void vector_resize(vector *v, size_t capacity)
{
    #ifdef DEBUG_ON
		printf("vector_resize: %d to %d\n", v->capacity, capacity);
    #endif

    int **data = realloc(v->data, sizeof(int *) * capacity);
	int *length = realloc(v->data, sizeof(int) * capacity);
	
    if (data)
    {
        v->data = data;
        v->capacity = capacity;
		v->length = length;
    }
}

void vector_add(vector *v, int *data, int length)
{
    if (v->capacity == v->size)
        vector_resize(v, v->capacity * 2);

    v->data[v->size] = malloc(length * sizeof(int));
    memcpy(v->data[v->size], data, length * sizeof(int));
    v->length[v->size] = length;
    v->size++;
}

void vector_set(vector *v, int index, int *data, int length)
{
    if (index >= 0 && index < v->size)
	{
        v->data[index] = data;
		v->length[index] = length;
	}
}

int *vector_get(vector *v, int index)
{
    if (index >= 0 && index < v->size)
        {
        	return v->data[index];
        }
    return NULL;
}

void vector_delete(vector *v, int index)
{
    if (index < 0 || index >= v->size)
        return;

    v->data[index] = NULL;
	int i; 
    for (i= 0; i < v->size - 1; i++)
    {
        v->data[i] = v->data[i + 1];
        v->data[i + 1] = NULL;
    }

    v->size--;

    if (v->size > 0 && v->size == v->capacity / 4)
        vector_resize(v, v->capacity / 2);
}

void vector_free(vector *v)
{
    free(v->data);
	v->data = NULL;
	v->capacity = v->size = 0;
}
