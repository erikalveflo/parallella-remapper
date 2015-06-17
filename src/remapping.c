#include "remapping.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

// A function pointer to a comparison function used when sorting a list.
typedef e_bool_t (*w_compare_func_t)(w_list_element_t, w_list_element_t);

////////////////////////////////////////////////////////////////////////////////

// Frees the data pointer to by [p] and sets to to NULL. Protects against
// multiple frees of the same data.
static void w_free(void **p)
{
	if (p != NULL && *p != NULL) {
		free(*p);
		*p = NULL;
	}
}

// Generates a random double between 0 and 1 using the supplied PRNG-state.
static double w_random01(pcg32_random_t *pcg)
{
	return (double)pcg32_random_r(pcg) / UINT32_MAX;
}

// Generates a random integer between 0 and [max] using the supplied PRNG-state.
static w_size_t w_random_range(w_size_t max, pcg32_random_t *pcg)
{
	return max * w_random01(pcg);
}

////////////////////////////////////////////////////////////////////////////////

w_coord_t w_get_coord(w_core_id_t id, w_size_t cols)
{
	w_coord_t coord;
	coord.row = id / cols;
	coord.col = id % cols;
	return coord;
}

w_core_id_t w_get_core_id(w_coord_t coord, w_size_t cols)
{
	return coord.col + coord.row * cols;
}

// Compute the Manhattan distance between two coordinates [a] and [b].
static double w_manhattan_distance(w_coord_t a, w_coord_t b)
{
	return abs(a.col - b.col) + abs(a.row - b.row);
}

////////////////////////////////////////////////////////////////////////////////

unsigned w_get_global_address(e_epiphany_t *device, w_core_id_t id, unsigned offset)
{
	w_coord_t coord = w_get_coord(id, device->cols);
	coord.col += device->base_coreid;
	return offset + coord.col * 0x00100000 + coord.row * 0x04000000;
}

w_error_code_t w_write(e_epiphany_t *device, w_core_id_t id, unsigned address, void *source, size_t size)
{
	w_coord_t coord = w_get_coord(id, device->cols);
	return e_write(device, coord.row, coord.col, address, source, size);
}

w_error_code_t w_read(e_epiphany_t *device, w_core_id_t id, unsigned address, void *target, size_t size)
{
	w_coord_t coord = w_get_coord(id, device->cols);
	return e_read(device, coord.row, coord.col, address, target, size);
}

////////////////////////////////////////////////////////////////////////////////

void w_init_list(w_list_t *list, w_size_t capacity)
{
	list->capacity = capacity;
	list->size = 0;

	if (capacity > 0) {
		list->elements = malloc(capacity * sizeof(list->elements[0]));
	} else {
		list->elements = NULL;
	}
}

void w_free_list(w_list_t *list)
{
	w_free((void **)&list->elements);

	list->capacity = 0;
	list->size = 0;
}

// Adds [element] to the end of [list].
static void w_push_back(w_list_t *list, w_list_element_t element)
{
	assert(list->capacity > list->size);

	list->elements[list->size] = element;
	++list->size;
}

// Removes element at [index] from [list].
static void w_remove_at(w_list_t *list, w_size_t index)
{
	assert(index < list->size);
	assert(list->size > 0);

	list->elements[index] = list->elements[list->size - 1];
	--list->size;
}

// Clears the list and grows its storage capacity.
static void w_clear_list(w_list_t *list, w_size_t new_capacity)
{
	if (new_capacity > list->capacity) {
		// Allocate more memory for the array.
		w_free_list(list);
		w_init_list(list, new_capacity);
	}

	list->size = 0;
}

// Copies the elements in [src] to [dst]. The lists can be of different sizes
// and capacities.
static void w_copy_list(w_list_t *dst, w_list_t *src)
{
	w_clear_list(dst, src->capacity);
	dst->size = src->size;
	memcpy(dst->elements, src->elements, src->capacity * sizeof(src->elements[0]));
}

// Sorts [list] based on the predicate [compare] where [compare] is a function
// which returns true when its first argument should go before its second
// argument. Implemented as bubble sort which is O(N^2).
static void w_sort(w_list_t *list, w_compare_func_t compare)
{
	w_size_t i, j;
	w_list_element_t tmp;

	for (i = 0; i < list->size; ++i) {
		for (j = 0; j < list->size; ++j) {
			if (compare(list->elements[i], list->elements[j]) == E_TRUE) {
				tmp = list->elements[i];
				list->elements[i] = list->elements[j];
				list->elements[j] = tmp;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////

w_matrix_t w_create_matrix(w_size_t rows, w_size_t cols)
{
	w_matrix_t matrix;
	uint16_t size = rows * cols;

	matrix.size = size;
	matrix.rows = rows;
	matrix.cols = cols;
	matrix.elements = calloc(size, sizeof(matrix.elements[0]));

	return matrix;
}

void w_free_matrix(w_matrix_t *matrix)
{
	w_free((void **)&matrix->elements);

	matrix->size = 0;
	matrix->rows = 0;
	matrix->cols = 0;
}

// Copies the elements in [src] to [dst]. The matrices must be of equal size.
static void w_copy_matrix(w_matrix_t *dst, w_matrix_t *src)
{
	assert(dst->size == src->size);

	memcpy(dst->elements, src->elements, dst->size * sizeof(dst->elements[0]));
}

void w_find_not_in_matrix(w_matrix_t *matrix, w_list_t *result, w_matrix_element_t needle)
{
	w_size_t i;

	w_clear_list(result, matrix->size);

	for (i = 0; i < matrix->size; ++i) {
		if (matrix->elements[i] != needle) {
			w_push_back(result, i);
		}
	}
}

void w_find_in_matrix(w_matrix_t *matrix, w_list_t *result, w_matrix_element_t needle)
{
	w_size_t i;

	w_clear_list(result, matrix->size);

	for (i = 0; i < matrix->size; ++i) {
		if (matrix->elements[i] == needle) {
			w_push_back(result, i);
		}
	}
}

void w_find_mask_in_matrix(w_matrix_t *matrix, w_list_t *result, w_matrix_element_t mask)
{
	w_size_t i;

	w_clear_list(result, matrix->size);

	for (i = 0; i < matrix->size; ++i) {
		if ((matrix->elements[i] & mask) != 0) {
			w_push_back(result, i);
		}
	}
}

void w_print_matrix(w_matrix_t *matrix, const char *name)
{
	w_size_t row, col;
	w_matrix_element_t e;

	printf("%s Matrix (rows: %u, cols: %u)\n", name, matrix->rows, matrix->cols);
	for (row = 0; row < matrix->rows; ++row) {
		printf("  ");
		for (col = 0; col < matrix->cols; ++col) {
			e = matrix->elements[col + row * matrix->cols];
			e != 0 ? (e == W_SIZE_MAX ? printf(" X ") : printf("%2u ", e)) : printf(" . ");
		}
		printf("\n");
	}
}

// Swaps element at index [a] with element at index [b] in [matrix].
static void w_swap_matrix_elements(w_matrix_t *matrix, uint16_t a, uint16_t b)
{
	assert(matrix->size > a && "index out of bounds");
	assert(matrix->size > b && "index out of bounds");

	w_matrix_element_t tmp = matrix->elements[a];
	matrix->elements[a] = matrix->elements[b];
	matrix->elements[b] = tmp;
}

// Shuffles the elements in [matrix] randomly using the suppled PRNG-state.
static void w_shuffle_matrix(w_matrix_t *matrix, pcg32_random_t *pcg)
{
	w_size_t i, a, b;
	w_size_t size;

	size = matrix->size;

	for (i = 0; i < 2 * size; ++i) {
		a = w_random_range(size, pcg);
		b = w_random_range(size, pcg);
		w_swap_matrix_elements(matrix, a, b);
	}
}

// Shuffles the [matrix] elements with [indices] randomly using the suppled
// PRNG-state.
static void w_shuffle_partial_matrix(w_matrix_t *matrix, w_list_t *indices, pcg32_random_t *pcg)
{
	w_size_t i, a, b;

	for (i = 2 * matrix->size; i > 0; --i) {
		a = w_random_range(indices->size, pcg);
		b = w_random_range(indices->size, pcg);

		a = indices->elements[a];
		b = indices->elements[b];

		w_swap_matrix_elements(matrix, a, b);
	}
}
////////////////////////////////////////////////////////////////////////////////

// Resets the mapping matrix to its original values: core 0 and location 0, core
// 1 at 1 ... core N at N.
static void w_reset_mapping_matrix(w_mapper_t *mapper)
{
	w_size_t i;

	for (i = 0; i < mapper->num_cores; ++i) {
		mapper->mapping_matrix.elements[i] = i;
	}
}

// Resets the constraint matrix to allow any task on all cores.
static void w_reset_constraint_matrix(w_mapper_t *mapper)
{
	w_size_t i;

	for (i = 0; i < mapper->num_cores; ++i) {
		mapper->constraint_matrix.elements[i] = W_ALLOW_ANY;
	}
}

w_mapper_t w_create_mapper(w_size_t rows, w_size_t cols)
{
	assert(rows > 0);
	assert(cols > 0);

	size_t size;
	uint16_t num_cores;
	w_mapper_t mapper;

	memset(&mapper, 0, sizeof(mapper));

	num_cores = rows * cols;
	mapper.num_cores = num_cores;
	mapper.rows = rows;
	mapper.cols = cols;

	size = num_cores * sizeof(mapper.links.indices[0]);
	mapper.links.indices = malloc(size);
	memset(mapper.links.indices, W_NO_LINK, size);

	mapper.allocation_matrix = w_create_matrix(rows, cols);
	mapper.constraint_matrix = w_create_matrix(rows, cols);
	mapper.fault_matrix      = w_create_matrix(rows, cols);
	mapper.assignment_matrix = w_create_matrix(rows, cols);
	mapper.mapping_matrix    = w_create_matrix(rows, cols);

	w_reset_constraint_matrix(&mapper);
	w_reset_mapping_matrix(&mapper);

	return mapper;
}

void w_free_mapper(w_mapper_t *mapper)
{
	mapper->num_cores = 0;

	mapper->links.capacity = 0;
	mapper->links.size = 0;
	w_free((void **)&mapper->links.indices);
	w_free((void **)&mapper->links.data);

	w_free_matrix(&mapper->allocation_matrix);
	w_free_matrix(&mapper->constraint_matrix);
	w_free_matrix(&mapper->fault_matrix);
	w_free_matrix(&mapper->assignment_matrix);
	w_free_matrix(&mapper->mapping_matrix);
}

void w_set_allocation_matrix(w_mapper_t *mapper, w_matrix_element_t *allocation_matrix)
{
	memcpy(mapper->allocation_matrix.elements, allocation_matrix,
		mapper->allocation_matrix.size * sizeof(w_matrix_element_t));
}

void w_set_constraint_matrix(w_mapper_t *mapper, w_matrix_element_t *constraint_matrix)
{
	memcpy(mapper->constraint_matrix.elements, constraint_matrix,
		mapper->constraint_matrix.size * sizeof(w_matrix_element_t));
}

void w_set_fault_matrix(w_mapper_t *mapper, w_matrix_element_t *fault_matrix)
{
	memcpy(mapper->fault_matrix.elements, fault_matrix,
		mapper->fault_matrix.size * sizeof(w_matrix_element_t));
}

void w_set_assignment_matrix(w_mapper_t *mapper, w_matrix_element_t *assignment_matrix)
{
	memcpy(mapper->assignment_matrix.elements, assignment_matrix,
		mapper->assignment_matrix.size * sizeof(w_matrix_element_t));
}

void w_add_link(w_mapper_t *mapper, w_core_id_t begin, w_core_id_t end)
{
	assert(mapper->num_cores > begin);
	assert(mapper->num_cores > end);

	w_link_collection_t *links;
	w_link_t link;
	w_size_t new_capacity;
	w_size_t insert_at, next;

	links = &mapper->links;

	if (links->size + 1 > links->capacity) {
		// Increase capacity to make room for more links.
		new_capacity = (links->capacity * 3) / 2;
		new_capacity = new_capacity < 8 ? 8 : new_capacity;

		printf("Increasing link storage capacity to %u\n", new_capacity);

		links->data = realloc(links->data, new_capacity * sizeof(links->data[0]));
		links->capacity = new_capacity;
	}

	// Insert at the end of the list. Use the next field to point to the
	// previous head of the array.
	insert_at = links->size;
	next = links->indices[begin];

	printf("Inserting link from %u to %u (at: %u, next: %u)\n", begin, end, insert_at, next);

	link.begin = begin;
	link.end = end;
	link.next = next;

	links->indices[begin] = insert_at;
	links->data[insert_at] = link;
	++links->size;
}

////////////////////////////////////////////////////////////////////////////////

// Check that the given [allocation]/assignment matrix upholds the constraints
// defined in the [constraint] matrix.
static w_error_code_t w_verify_constraints(w_matrix_t *constraint, w_matrix_t *allocation)
{
	assert(constraint->size == allocation->size);

	w_size_t i;
	w_task_id_t task;

	for (i = 0; i < allocation->size; ++i) {
		task = allocation->elements[i];
		if (task != 0 && (constraint->elements[i] & task) == 0) {
			printf("Constraint violation: task %u not allowed on core %u\n", task, i);
			return E_ERR;
		}
	}

	return E_OK;
}

// Check that the given assignment matrix does not use faulty cores from the
// [faults] matrix.
static w_error_code_t w_verify_assignment(w_matrix_t *faults, w_matrix_t *assignment)
{
	assert(faults->size == assignment->size);

	w_size_t i;
	w_task_id_t task;

	for (i = 0; i < assignment->size; ++i) {
		task = assignment->elements[i];
		if (task != 0 && faults->elements[i] != 0) {
			printf("Faulty core %u assigned to task %u\n", i, task);
			return E_ERR;
		}
	}

	return E_OK;
}

// Counts the number of set bits (ones) in the supplied 32-bit number. Adopted
// from code by Matt Howells http://stackoverflow.com/a/109025/1009118
static uint32_t w_count_high_bits(uint32_t i)
{
     i = i - ((i >> 1) & 0x55555555);
     i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
     return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

// A predicate for sorting by the number of set bits in ascending order.
static e_bool_t w_sort_by_number_of_bits_asc(w_list_element_t a, w_list_element_t b)
{
	return w_count_high_bits(a) > w_count_high_bits(b) ? E_TRUE : E_FALSE;
}

w_error_code_t w_assign_naive(w_mapper_t *mapper)
{
	w_error_code_t error_code;

	w_size_t i, j;

	w_list_t faults, spares, choices;
	w_core_id_t faulty, spare;
	w_task_id_t task;
	e_bool_t spare_found;
	e_bool_t available_cores[mapper->num_cores];

	w_init_list(&faults, mapper->num_cores);
	w_init_list(&spares, mapper->num_cores);
	w_init_list(&choices, mapper->num_cores);

	if (w_verify_constraints(&mapper->constraint_matrix, &mapper->allocation_matrix) != E_OK) {
		printf("WARNING: Constraint violation in allocation matrix.\n");
	}

	// Reset mapper: clear changes and use the allocation matrix as the initial
	// solution.
	w_reset_mapping_matrix(mapper);
	w_copy_matrix(&mapper->assignment_matrix, &mapper->allocation_matrix);

	w_find_not_in_matrix(&mapper->fault_matrix, &faults, 0);
	if (faults.size == 0) {
		// No faults: the allocation matrix is a good assignment matrix.
		goto success;
	}

	// Faults are present in the system. Find the spares.
	w_find_in_matrix(&mapper->allocation_matrix, &spares, 0);

	if (faults.size > spares.size) {
		printf("ERROR: Too many faults. Only %u spare cores are available but %u are faulty.\n", spares.size, faults.size);
		goto error;
	}

	// Mark all cores as either available or not. Allocated and faulty cores are
	// not available.
	for (i = 0; i < mapper->num_cores; ++i) {
		if (mapper->fault_matrix.elements[i] == 0 && mapper->allocation_matrix.elements[i] == 0) {
			available_cores[i] = E_TRUE;
		} else {
			available_cores[i] = E_FALSE;
		}
	}

	// This is the naive part. Replace each faulty core with the first spare
	// core.
	for (i = 0; i < faults.size; ++i) {
		faulty = faults.elements[i];
		task = mapper->allocation_matrix.elements[faulty];

		if (task == 0)
			continue;

		// Find the cores form which a spare can be chosen. The choices are
		// limited by constraints and faults. It should be possible to cache
		// these results.
		w_find_mask_in_matrix(&mapper->constraint_matrix, &choices, task);

		// Sort the cores in order of restrictiveness. This will assign
		// restrictive cores first; for example, if core 1 only accepts task 1
		// then it will be assigned prior to core 0 which accepts any task.
		w_sort(&choices, &w_sort_by_number_of_bits_asc);

		// Find an unassigned core among the choices that fulfill the
		// constraints.
		spare_found = E_FALSE;
		for (j = 0; j < choices.size; ++j) {
			spare = choices.elements[j];
			if (available_cores[spare] == E_TRUE) {
				spare_found = E_TRUE;
				break;
			}
		}

		if (spare_found == E_FALSE) {
			printf("ERROR: Spare not found for core %u. Task %u is constrained to %u location(s) of %u free spare location(s).\n", faulty, task, choices.size, spares.size);
			w_print_matrix(&mapper->assignment_matrix, "Partial Assignment");
			goto error;
		}

		printf("Fault at %u replaced by %u\n", faulty, spare);
		mapper->assignment_matrix.elements[spare] = task;
		mapper->assignment_matrix.elements[faulty] = 0;
		available_cores[spare] = E_FALSE;

		// Keep track of all changes.
		w_swap_matrix_elements(&mapper->mapping_matrix, faulty, spare);
	}

	if (w_verify_constraints(&mapper->constraint_matrix, &mapper->assignment_matrix) != E_OK) {
		printf("ERROR: Constraint violation in assignment matrix.\n");
		goto error;
	}

	// The assignment was a success.
	goto success;

success:
	error_code = E_OK;
	goto cleanup;

error:
	error_code = E_ERR;
	goto cleanup;

cleanup:
	w_free_list(&faults);
	w_free_list(&spares);
	w_free_list(&choices);
	return error_code;
}

////////////////////////////////////////////////////////////////////////////////

// Computes the probability of accepting a neighbor solution when running
// simulated annealing for the given temperature and energies. The probability
// is a number between 0 and 1; where 1 means that the neighbor solution should
// be accepted.
static double w_acceptance_propability(double current_energy, double neighbor_energy, double temperature)
{
	if (neighbor_energy < current_energy)
		return 1.0;

	return exp((current_energy - neighbor_energy) / temperature);
}

// Computes the energy of the given [solution]. The energy is the cost function
// that simulated annealing tries to minimize. A high energy equals a high cost
// and therefore a bad solution.
//
// The position of each core is given by [solution]. solution[0] is the position
// of the zeroth core. If, for example, solution[0] == 2 core zero can be found
// at position 2.
static double w_compute_energy(w_mapper_t *mapper, w_matrix_t *solution)
{
	double energy;
	w_link_t *link;
	w_size_t i;
	w_coord_t begin, end;

	energy = 0.0;

	// Compute the distance between each connected core in the current solution.
	for (i = 0; i < mapper->links.size; ++i) {
		link = &mapper->links.data[i];

		begin = w_get_coord(solution->elements[link->begin], mapper->cols);
		end = w_get_coord(solution->elements[link->end], mapper->cols);

		energy += w_manhattan_distance(begin, end);
	}

	return energy;
}

// Returns true if the cores [a] and [b] are allowed to swap positions based on
// the tasks allocated to those cores and their constraints.
static e_bool_t w_allow_swap(w_mapper_t *mapper, w_matrix_t *solution, w_core_id_t a, w_core_id_t b)
{
	w_core_id_t pos_a, pos_b;
	w_task_id_t task_a, task_b;
	w_constraint_mask_t mask_a, mask_b;
	w_matrix_element_t fault_a, fault_b;

	pos_a = solution->elements[a];
	pos_b = solution->elements[b];

	task_a = mapper->allocation_matrix.elements[a];
	task_b = mapper->allocation_matrix.elements[b];

	mask_a = mapper->constraint_matrix.elements[pos_a];
	mask_b = mapper->constraint_matrix.elements[pos_b];

	fault_a = mapper->fault_matrix.elements[pos_a];
	fault_b = mapper->fault_matrix.elements[pos_b];

	if (task_a != 0) {
		if (fault_b != 0 || (mask_b & task_a) == 0) {
			return E_FALSE;
		}
	}

	if (task_b != 0) {
		if (fault_a != 0 || (mask_a & task_b) == 0) {
			return E_FALSE;
		}
	}

	return E_TRUE;
}

w_sa_config_t w_create_sa_config()
{
	w_sa_config_t c;
	c.start_temperature = 1.0e3;
	c.stop_temperature = 1.0e-10;
	c.cooling_rate = 0.3e-3;
	c.partial = E_FALSE;
	c.seed = 0;
	return c;
}

// Get two indices [pa] and [pb] to swap. These indices are based on [choices].
// Returns true if the returned indices represent an allowed swap (that does not
// violate design constraints.)
static e_bool_t w_get_swap_indices(
	w_mapper_t *mapper,
	w_list_t *choices,
	w_list_t *swappable,
	w_matrix_t *current_solution,
	w_matrix_t *inverse_solution,
	pcg32_random_t *pcg,
	w_size_t *pa, // Out argument. Index in [current_solution].
	w_size_t *pb  // Out argument. Index in [current_solution].
) {
	w_size_t i, a, b;
	w_core_id_t id;

	// Create a list of swappable indices. [choices] contain core IDs but we
	// need the position of the cores (their index in [current_solution].)
	swappable->size = choices->size;
	for (i = 0; i < choices->size; ++i) {
		id = choices->elements[i];
		swappable->elements[i] = inverse_solution->elements[id];
	}

	a = swappable->elements[w_random_range(swappable->size, pcg)];
	b = swappable->elements[w_random_range(swappable->size, pcg)];

	*pa = a;
	*pb = b;

	// Does the swap violate constraints?
	return w_allow_swap(mapper, current_solution, a, b);
}

// Swap index [a] and [b] in [current_solution]. Call this function again (with
// the same arguments) to undo the swap.
static void w_swap(w_matrix_t *current_solution, w_matrix_t *inverse_solution, w_size_t a, w_size_t b)
{
	w_swap_matrix_elements(current_solution, a, b);
	w_swap_matrix_elements(inverse_solution, current_solution->elements[a], current_solution->elements[b]);
}

w_error_code_t w_assign_sa(w_mapper_t *mapper, w_sa_config_t *_config)
{
	w_sa_config_t config = _config == NULL ? w_create_sa_config() : *_config;

	w_error_code_t error_code;

	double temperature = config.start_temperature;
	const double cooling_factor = 1.0 - config.cooling_rate;

	// Expected number of iterations under ideal conditions.
	const uint32_t expected_iterations = ceil((log(config.stop_temperature) - log(config.start_temperature)) / log(cooling_factor));

	// If the number of failed iterations exceed this number, no solution was
	// found in a reasonable amount of time.
	const uint32_t max_failed_iterations = expected_iterations;

	uint32_t iteration_count = 0; // Number of performed iterations.
	uint32_t failure_count = 0; // Number of iterations that violated constraints.

	w_size_t i, a, b, position, task; // Temporary variables.
	w_core_id_t id;

	double accept = 0.0; // Acceptance probability.

	double neighbor_energy = 0.0; // Energy of new solution, not yet accepted.
	double current_energy = 0.0; // Energy of currently accepted solution.
	double lowest_energy = 0.0; // Lowest found energy (see lowest_solution.)

	e_bool_t allow_swap; // Is the swap allowed? Does it violate constraints?

	// These matrices holds the position of each core (same as
	// [mapper.mapping_matrix].) For example, element 0 holds the position of
	// core 0. Translates core ID to grid position.
	w_matrix_t current_solution = w_create_matrix(mapper->rows, mapper->cols);
	w_matrix_t lowest_solution = w_create_matrix(mapper->rows, mapper->cols);

	// This matrix holds the core ID at each position. For example, element 0
	// holds the core ID at position 0. Translates grid position to core ID.
	w_matrix_t inverse_solution = w_create_matrix(mapper->rows, mapper->cols);

	// Initialize the random number generator. We could use time() here but
	// we'll let the user decide if that is useful.
	pcg32_random_t pcg = PCG32_INITIALIZER;
	pcg.state += config.seed + 2;

	// Lists all core IDs that we are allowed to swap. For partial remapping
	// this list is small; otherwise it contains all core IDs.
	w_list_t choices;
	w_init_list(&choices, mapper->num_cores);

	// Lists all indices in [current_solution] that we are allowed to swap. This
	// list is based on [choices] but changes with each swap.
	w_list_t swappable;
	w_init_list(&swappable, mapper->num_cores);

	// Temporary list of all spare cores in the assignment matrix.
	w_list_t spares;
	w_init_list(&spares, 0);

	if (w_verify_constraints(&mapper->constraint_matrix, &mapper->allocation_matrix) != E_OK) {
		printf("WARNING: Constraint violation in allocation matrix.\n");
	}

	if (config.partial == E_FALSE) {
		// All cores are assigned during non-partial assignment.
		for (i = 0; i < mapper->num_cores; ++i) {
			w_push_back(&choices, i);
		}

		// Reset and forget about previous solutions.
		w_reset_mapping_matrix(mapper);
		w_copy_matrix(&mapper->assignment_matrix, &mapper->allocation_matrix);
	} else {
		// Only a limited number of cores are assigned during partial
		// assignment. These cores are added to [choices]. We only need to remap
		// assigned faulty cores.

		// TODO: Test if a valid assignment matrix exists.
		// NOTE: [spares] is a temporary list only used here.

		// Add all faulty cores to the list of choices, and then remove all
		// faulty non-assigned cores.
		w_find_not_in_matrix(&mapper->fault_matrix, &choices, 0);
		for (i = 0; i < choices.size; ++i) {
			id = choices.elements[i];
			if (mapper->assignment_matrix.elements[id] == 0) {
				w_remove_at(&choices, i);
				--i;
			}
		}

		// Add all spare cores to the list of choices.
		w_find_in_matrix(&mapper->assignment_matrix, &spares, 0);
		for (i = 0; i < spares.size; ++i) {
			id = spares.elements[i];
			if (mapper->fault_matrix.elements[id] == 0) {
				w_push_back(&choices, id);
			}
		}
	}

	// Operate on [current_solution] instead of [mapper.mapping_matrix] even
	// though they represent the same thing.
	w_copy_matrix(&current_solution, &mapper->mapping_matrix);

	// Construct the inverse of the [current_solution].
	for (i = 0; i < mapper->num_cores; ++i) {
		position = current_solution.elements[i];
		inverse_solution.elements[position] = i;
	}

	// Randomize the current solution.
	for (i = 2 * mapper->num_cores; i > 0; ++i) {
		w_get_swap_indices(mapper, &choices, &swappable,
			&current_solution, &inverse_solution, &pcg, &a, &b);
		w_swap(&current_solution, &inverse_solution, a, b);
	}

	// Store best known lowest-energy solution, its energy and the energy of the
	// current solution.
	w_copy_matrix(&lowest_solution, &current_solution);
	current_energy = w_compute_energy(mapper, &current_solution);
	lowest_energy = current_energy;

	printf("Start energy %.2f (expected iterations: %u)\n", lowest_energy, expected_iterations);

	// Simulate annealing until the target temperature is reached.
	while (temperature > config.stop_temperature) {
		++iteration_count;

		// Swap the position of two random cores and test if that neighboring
		// solution is a better solution than the current one.
		allow_swap = w_get_swap_indices(mapper, &choices, &swappable,
			&current_solution, &inverse_solution, &pcg, &a, &b);

		// Does the swap violate constraints?
		if (allow_swap == E_FALSE) {
			// Abort if no solution can be found within a reasonable amount of
			// time.
			++failure_count;
			if (failure_count > max_failed_iterations) {
				printf("ERROR: Unable to find solution after %u failed iterations.\n", failure_count);
				goto error;
			}
			continue;
		}

		// Perform the swap and compute the new energy.
		w_swap(&current_solution, &inverse_solution, a, b);
		neighbor_energy = w_compute_energy(mapper, &current_solution);

		// If the neighbor solution is worse than the current one there is still
		// a chance that it will be accepted. Compute that chance.
		accept = w_acceptance_propability(current_energy, neighbor_energy,
			temperature);

		if (accept < w_random01(&pcg)) {
			// Don't accept the solution.
			w_swap(&current_solution, &inverse_solution, a, b);
		} else {
			// Accept the solution.
			current_energy = neighbor_energy;

			// Did we find a lower energy solution than previously recorded?
			if (lowest_energy > current_energy) {
				lowest_energy = current_energy;
				w_copy_matrix(&lowest_solution, &current_solution);
				printf("  energy: %.2f, iterations: %u, temperature: %.2f\n", current_energy, iteration_count, temperature);
			}
		}

		// Reduce the temperature slowly.
		temperature *= cooling_factor;
	}

	printf("End energy %.2f (iterations: %u, failed: %u)\n", lowest_energy, iteration_count, failure_count);

	// Copy the best solution to a publicly accessible matrix.
	w_copy_matrix(&mapper->mapping_matrix, &lowest_solution);

	// Apply the solution.
	for (i = 0; i < mapper->num_cores; ++i) {
		position = mapper->mapping_matrix.elements[i];
		task = mapper->allocation_matrix.elements[i];
		mapper->assignment_matrix.elements[position] = task;
	}

	if (w_verify_constraints(&mapper->constraint_matrix, &mapper->assignment_matrix) != E_OK) {
		printf("ERROR: Constraint violation in assignment matrix.\n");
		goto error;
	}

	if (w_verify_assignment(&mapper->fault_matrix, &mapper->assignment_matrix) != E_OK) {
		printf("ERROR: Faulty cores assigned in assignment matrix.\n");
		goto error;
	}

	// The assignment was a success.
	goto success;

success:
	error_code = E_OK;
	goto cleanup;

error:
	error_code = E_ERR;
	goto cleanup;

cleanup:
	w_free_list(&choices);
	w_free_list(&swappable);
	w_free_list(&spares);
	return error_code;
}

void w_set_partial_constraint_pattern(w_mapper_t *mapper)
{
	const int num_spare_cols = 1;

	int row, col;
	w_coord_t coord;
	w_core_id_t id;

	for (col = num_spare_cols; col < mapper->cols; col += num_spare_cols + 1) {
		coord.col = col;
		for (row = 0; row < mapper->rows; ++row) {
			coord.row = row;
			id = w_get_core_id(coord, mapper->cols);
			mapper->constraint_matrix.elements[id] = 0;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////

w_error_code_t w_load(e_epiphany_t *device, w_list_t *core_ids, char *executable)
{
	const e_bool_t start_program = E_FALSE;

	int i;
	w_error_code_t status;
	w_coord_t coord;

	for (i = 0; i < core_ids->size; ++i) {
		coord = w_get_coord(core_ids->elements[i], device->cols);
		status = e_load(executable, device, coord.row, coord.col, start_program);
		if (status != E_OK)
			return status;
	}

	return E_OK;
}
