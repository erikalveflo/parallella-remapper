#pragma once

#include "pcg.h"

#include <e-hal.h>
#include <e-loader.h>

// Used for documentation purposes to brand functions that can fail.
typedef int w_error_code_t;
typedef uint16_t w_matrix_element_t;
typedef uint16_t w_list_element_t;
typedef uint16_t w_core_id_t;
typedef uint16_t w_task_id_t;
typedef uint16_t w_constraint_mask_t;
typedef uint16_t w_size_t;

enum {
	W_SIZE_MAX = (w_size_t)~0,

	// Used in contraint matrices to allow any task (unconstrained core).
	W_ALLOW_ANY = (w_matrix_element_t)~0,

	// Indicates that a [w_mapper_t.link_index] is unused.
	W_NO_LINK = W_SIZE_MAX,
};

////////////////////////////////////////////////////////////////////////////////

// A w_matrix_t is a matrix stored as a one-dimensional array. Always free with
// w_free_matrix() when done.
typedef struct {
	w_size_t size; // Number of elements.
	w_size_t rows; // Number of rows.
	w_size_t cols; // Number of columns.
	w_matrix_element_t *elements; // Array of matrix elements.
} w_matrix_t;

// A w_list_t is a list/array with known capacity and size. Contrary to ordinary
// C-arrays this list keeps track of its capacity and the number of elements
// currently occupied. Always initialize with w_init_list() prior to use and
// always free with w_free_list() when done.
typedef struct {
	w_size_t capacity; // Number of elements the list can hold.
	w_size_t size; // Number of elements currently occupied.
	w_list_element_t *elements; // Array of list elements.
} w_list_t;

// A w_link_t represents a data link/connection between two cores. Connected
// cores send data to each other. Links are stored as a linked list where [next]
// is the index of the next element in the linked list.
typedef struct {
	w_core_id_t begin; // Link begins at this core.
	w_core_id_t end; // Link ends at this core.
	w_size_t next; // Index of the next item in the linked list or W_NO_LINK.
} w_link_t;

// A w_link_collection_t is a dynamic collection of w_link_t. Links are stored
// as a linked list where [indices] points to the first element in the list and
// [data] holds all links in the collection. [indices] includes an element for
// each core which defaults to W_NO_LINK for cores without links.
typedef struct {
	w_size_t capacity; // Number of links [data] can hold.
	w_size_t size; // Number of occupied links in [data].
	w_size_t *indices; // Maps a core ID to an element in [data] or W_NO_LINK.
	w_link_t *data; // A list of links.
} w_link_collection_t;

// All state information needed by the mapping system.
typedef struct {
	w_size_t num_cores; // Number of cores in the device.
	w_size_t rows; // Number of rows in the device.
	w_size_t cols; // Number of columns in the device.

	w_link_collection_t links; // A collection of data links.

	// The allocation matrix defines the original (intended) location of all
	// tasks and the number of cores allocated for each task. For example, in
	// matrix [1 1; 0 0], a task 1 is allocated to core 0 and 1. Tasks should
	// use power of two numbers, such as 1, 2, 4, 8, ..., 2^N. Zero denotes
	// spare cores.
	w_matrix_t allocation_matrix;

	// The constraints matrix is used to constrain tasks to specific cores.
	// Tasks are either allowed on all cores or constrained to a specific number
	// of cores defined in the constraint matrix. The constraints are specified
	// as a bit mask. The matrix [0 0; 1 0] constrain task 1 to core 3; it does
	// not allow other tasks.
	w_matrix_t constraint_matrix;

	// Non-zero elements indicates faulty cores. For example, in matrix [1 0; 0
	// 0] the zeroth core is faulty.
	w_matrix_t fault_matrix;

	// The assignment matrix maps a task to a core. This is the result of the
	// assignment algorithm. For example, in matrix [1 0; 2 2], core 0 should
	// load task 1, and core 2 and 3 should load task 2.
	w_matrix_t assignment_matrix;

	// The mapping matrix indicates the position of each core. For example,
	// mapping_matrix[0] is the position of the zeroth core.
	w_matrix_t mapping_matrix;
} w_mapper_t;

// A two-dimensional coordinate.
typedef struct {
	w_core_id_t row; // Row.
	w_core_id_t col; // Column.
} w_coord_t;

// Configuration values for the simulated annealing assignment.
typedef struct {
	double stop_temperature; // Initial temperature.
	double start_temperature; // Temperature when to stop.
	double cooling_rate; // Temperature cooling rate.
	e_bool_t partial; // Perform partial remapp?

	// Seed for the pseudo random number generator. Use time() to generate a
	// unique solution each second.
	uint32_t seed;
} w_sa_config_t;

////////////////////////////////////////////////////////////////////////////////

// Convert a core ID to coordinates given the number of columns in the device.
w_coord_t w_get_coord(w_core_id_t id, w_size_t cols);

// Convert coordinates to a core ID given the number of columns in the device.
w_core_id_t w_get_core_id(w_coord_t coord, w_size_t cols);

////////////////////////////////////////////////////////////////////////////////

// Returns the global address for a core (row, col). Host processor
// implementation of e_get_global_address.
unsigned w_get_global_address(e_epiphany_t *device, w_core_id_t id, unsigned offset);

// Writes [size] number of bytes from [source] to device memory at [address] on
// core [id].
w_error_code_t w_write(e_epiphany_t *device, w_core_id_t id, unsigned address, void *source, size_t size);

// Reads [size] number of bytes to [target] from device memory at [address] on
// core [id].
w_error_code_t w_read(e_epiphany_t *device, w_core_id_t id, unsigned address, void *target, size_t size);

////////////////////////////////////////////////////////////////////////////////

// Initializes [list] with the capacity to hold the given number of elements.
// Always call this function prior to using the list and prior to calling
// w_free_list().
void w_init_list(w_list_t *list, w_size_t capacity);

// Frees the memory allocated for [list].
void w_free_list(w_list_t *list);

////////////////////////////////////////////////////////////////////////////////

// Returns a new matrix of size [rows] and [cols].
w_matrix_t w_create_matrix(w_size_t rows, w_size_t cols);

// Frees the memory allocated for [matrix].
void w_free_matrix(w_matrix_t *matrix);

// Fills [result] with the indices of all elements in [matrix] not equal to
// [needle].
void w_find_not_in_matrix(w_matrix_t *matrix, w_list_t *result, w_matrix_element_t needle);

// Fills [result] with the indices of all elements in [matrix] equal to
// [needle].
void w_find_in_matrix(w_matrix_t *matrix, w_list_t *result, w_matrix_element_t needle);

// Fills [result] with the indices of all elements in [matrix] with any bit in
// common with [mask].
void w_find_mask_in_matrix(w_matrix_t *matrix, w_list_t *result, w_matrix_element_t mask);

// Prints [matrix] to stdout.
void w_print_matrix(w_matrix_t *matrix, const char *name);

////////////////////////////////////////////////////////////////////////////////

// Creates a new mapper for a device with the supplied number of [rows] and
// [cols].
w_mapper_t w_create_mapper(w_size_t rows, w_size_t cols);

// Frees the memory allocated for the supplied [mapper].
void w_free_mapper(w_mapper_t *mapper);

// Sets the allocation matrix of [mapper]. The number of elements in the
// supplied matrix should correspond to the number of cores in the device. The
// matrix elements are copied into a matrix owned by [mapper].
void w_set_allocation_matrix(w_mapper_t *mapper, w_matrix_element_t *allocation_matrix);

// Sets the constraint matrix of [mapper]. The number of elements in the
// supplied matrix should correspond to the number of cores in the device. The
// matrix elements are copied into a matrix owned by [mapper].
void w_set_constraint_matrix(w_mapper_t *mapper, w_matrix_element_t *constraint_matrix);

// Sets the fault matrix of [mapper]. The number of elements in the
// supplied matrix should correspond to the number of cores in the device. The
// matrix elements are copied into a matrix owned by [mapper].
void w_set_fault_matrix(w_mapper_t *mapper, w_matrix_element_t *fault_matrix);

// Sets the assignment matrix of [mapper]. The number of elements in the
// supplied matrix should correspond to the number of cores in the device. The
// matrix elements are copied into a matrix owned by [mapper].
void w_set_assignment_matrix(w_mapper_t *mapper, w_matrix_element_t *assignment_matrix);

// Adds a data link/connection from [begin] core to [end] core. The [begin] core
// sends data to the [end] core. The link defines a connection between the two
// cores which is used during assignment.
void w_add_link(w_mapper_t *mapper, w_core_id_t begin, w_core_id_t end);

////////////////////////////////////////////////////////////////////////////////

// Assigns a task to each core using a naive approach. This approach is fast but
// does not consider links or other optimizations.
w_error_code_t w_assign_naive(w_mapper_t *mapper);

////////////////////////////////////////////////////////////////////////////////

// Returns the default simulated annealing configuration. These defaults will
// not work well for all applications and should be tweaked to suite your needs.
w_sa_config_t w_create_sa_config();

// Assigns a task to each core using simulated annealing. This approach is slow
// but tries to optimize network distances, network utilization, and other
// metrics defined by a cost function. This approach only works when links have
// been added using w_add_link().
//
// Configuration values are sent with the [config] pointer. If NULL the default
// values returned by w_create_sa_config() are used.
w_error_code_t w_assign_sa(w_mapper_t *mapper, w_sa_config_t *config);

// Update the constraint matrix to reserve columns of spare cores for use with
// partial remapping. Every second column is reserved, limiting the maximum
// distance to a spare core to one. The reserved columns are 2, 5, 8, ..., 3+2n.
void w_set_partial_constraint_pattern(w_mapper_t *mapper);

////////////////////////////////////////////////////////////////////////////////

// Loads the [executable] onto each core in the list of [core_ids].
w_error_code_t w_load(e_epiphany_t *device, w_list_t *core_ids, char *executable);
