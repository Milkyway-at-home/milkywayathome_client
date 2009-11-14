#include "simulation_heap.h"
#include "search_parameters.h"
#include "stdlib.h"
#include "limits.h"
#include "float.h"

double current_time;
int heap_size, max_heap_size;
HEAP_NODE **heap;

void simulation_heap__init(int max_size, int number_parameters, char *time_template_filename) {
	int i;

	current_time = 0;
	heap_size = 0;
	max_heap_size = max_size + 1;
	heap = (HEAP_NODE**)malloc(sizeof(HEAP_NODE*) * max_heap_size);

	for (i = 0; i < max_heap_size; i++) {
		heap[i] = (HEAP_NODE*)malloc(sizeof(HEAP_NODE));
		heap[i]->time = DBL_MAX;
		init_search_parameters(&(heap[i]->sp), number_parameters);
	}
	printf("created heap of size %d, using time template file: %s\n", max_size, time_template_filename);
}


int simulation_heap__parent(int position) {
	return (position - 1) / 2;
}

int simulation_heap__left_child(int position) {
	return (2 * position) + 1;
}

int simulation_heap__right_child(int position) {
	return (2 * position) + 2; 
}


double get_simulation_time() {
	return current_time + drand48();
}


void set_heap_node(int position, SEARCH_PARAMETERS *sp) {
	copy_search_parameters__no_alloc(sp, heap[position]->sp);
	heap[position]->time = get_simulation_time();
}

void remove_heap_node(int position, SEARCH_PARAMETERS *sp) {
	copy_search_parameters__no_alloc(heap[position]->sp, sp);
	heap[position]->time = DBL_MAX;
}

void swap_heap_node(int position1, int position2) {
	HEAP_NODE *temp;
	temp = heap[position1];
	heap[position1] = heap[position2];
	heap[position2] = temp;
}


void simulation_heap__insert(SEARCH_PARAMETERS *sp) {
	int current, parent;

	current = heap_size;
//	printf("inserting, current: %d\n", current);
	parent = simulation_heap__parent(current);

	set_heap_node(current, sp);

	while (current != 0 && heap[parent]->time > heap[current]->time) {
//		printf("swapping %d with %d\n", parent, current);
//		printf("times: %d with %d\n", heap[parent]->time, heap[current]->time);
		swap_heap_node(parent, current); 
		current = parent;
		parent = simulation_heap__parent(current);
	}
	heap_size++;
}

void simulation_heap__remove_min(SEARCH_PARAMETERS *sp) {
	int left, right, current;

	current_time = heap[0]->time;
	printf("[time: %lf] ", current_time);
//	printf("removing head node\n");
	remove_heap_node(0, sp);

//	printf("swapping head node\n");
	swap_heap_node(0, heap_size-1);

	current = 0;
	left = simulation_heap__left_child(current);
	right = simulation_heap__right_child(current);
//	printf("current: %d [%d], left: %d [%d], right: %d [%d]\n", current, heap[current]->time, left, heap[left]->time, right, heap[right]->time);

	while (heap[current]->time > heap[left]->time || heap[current]->time > heap[right]->time) {
//		printf("looping in remove\n");
		if (heap[left]->time > heap[right]->time) {
			swap_heap_node(current, right);
			current = right;
		} else {
			swap_heap_node(current, left);
			current = left;
		}
		left = simulation_heap__left_child(current);
		right = simulation_heap__right_child(current);

		if (left >= max_heap_size) {
			break;
		}
		if (right >= max_heap_size) {
			if (left < max_heap_size) {
				if (heap[current]->time > heap[left]->time) {
					printf("swapping current and left\n");
					swap_heap_node(current, left);
				}
			}
			break;
		} 
//		printf("current: %d [%d], left: %d [%d], right: %d [%d]\n", current, heap[current]->time, left, heap[left]->time, right, heap[right]->time);
	}
	heap_size--;
}
