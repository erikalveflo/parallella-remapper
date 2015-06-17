//
// Links example.
//
// Copyright (c) 2015, Erik Alveflo.
//

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <e-hal.h>
#include <e-loader.h>

#include <remapping.h>
#include "shared.h"

// Allocate multiple cores for task 1. Use power of two (2, 4, 8...) when adding
// additional tasks. A task is an Epiphany device program.
w_matrix_element_t allocation_matrix[16] = {
	1, 1, 1, 1, // Element 0 corresponds to core (0,0), element 1 to (0,1), and
	0, 0, 0, 0, // element 4 to (1,0) etc.
	0, 0, 0, 0,
	0, 0, 0, 0,
};

// Allow task 1 on all cores. This is a per core bit-mask.
#define X W_ALLOW_ANY
w_matrix_element_t constraint_matrix[16] = {
	X, X, X, X,
	X, X, X, X,
	X, X, X, X,
	X, X, X, X,
};
#undef X

// Simulate that a core is faulty by setting its element to 1. Task 1 will be
// assigned to a healthy core (any core with a zero.)
w_matrix_element_t fault_matrix[16] = {
	1, 0, 0, 0,
	0, 0, 0, 0,
	0, 0, 0, 0,
	0, 0, 0, 0,
};

static void msleep(int ms)
{
	usleep(1000 * ms);
}

static void connect_task1(w_mapper_t *mapper, w_list_t *tasks)
{
	int i;
	w_core_id_t core_id, next_id;

	for (i = 0; i < tasks->size; ++i) {
		core_id = tasks->elements[i];
		next_id = tasks->elements[(i + 1) % tasks->size];
		w_add_link(mapper, core_id, next_id);
	}
}

// Initializes and setup all cores assigned to execute task 1.
static void init_task1(e_epiphany_t *device, w_list_t *tasks)
{
	int i;
	w_core_id_t core_id, next_id;
	Mailbox mailbox;

	for (i = 0; i < tasks->size; ++i) {
		core_id = tasks->elements[i];
		next_id = tasks->elements[(i + 1) % tasks->size];

		// This core should have a pointer to the next core's mailbox.
		memset(&mailbox, 0, sizeof(mailbox));
		mailbox.next_mailbox = w_get_global_address(device, next_id, _MAILBOX_ADDRESS);

		w_write(device, core_id, _MAILBOX_ADDRESS, &mailbox, sizeof(mailbox));
	}
}

int main(int argc, char *argv[])
{
	e_platform_t platform;
	e_epiphany_t device;

	w_mapper_t mapper;
	w_sa_config_t sa_config;
	w_list_t task1;
	w_core_id_t first_id, last_id;

	Mailbox mailbox;
	memset(&mailbox, 0, sizeof(mailbox));

	w_init_list(&task1, 0);

	printf("=== Initializing system\n");
	e_set_host_verbosity(H_D0);

	e_init(NULL);
	e_reset_system();
	e_get_platform_info(&platform);

	printf("=== Creating workgroup\n");
	e_open(&device, 0, 0, platform.rows, platform.cols);

	printf("=== Mapping device program\n");
	e_reset_group(&device);

	// Initialize the mapping system: we will use this to automatically map our
	// application, to assign a device program to each core.
	mapper = w_create_mapper(platform.rows, platform.cols);

	// Tell the mapper about our application. It needs to know about allocated
	// tasks, constraints and faults are optional. By default the mapper assumes
	// no faults and no constraints.
	w_set_allocation_matrix(&mapper, allocation_matrix);
	w_set_constraint_matrix(&mapper, constraint_matrix);
	w_set_fault_matrix(&mapper, fault_matrix);

	// Find the ID of all cores allocated for task 1, and create a link between
	// each core. Links are used to indicate which tasks that communicate with
	// each other, and to minimize the distance between communicating tasks.
	w_find_in_matrix(&mapper.allocation_matrix, &task1, 1);
	connect_task1(&mapper, &task1);

	// Map the application using simulated annealing. This is will optimize our
	// poor allocation.
	sa_config = w_create_sa_config();
	if (w_assign_sa(&mapper, &sa_config) != E_OK) {
		printf("ERROR: Assignment failed.\n");
		return 1;
	}

	w_print_matrix(&mapper.assignment_matrix, "Assignment");
	w_print_matrix(&mapper.mapping_matrix, "Mapping");

	// Find the ID of all cores assigned to task 1.
	w_find_in_matrix(&mapper.assignment_matrix, &task1, 1);
	first_id = task1.elements[0];
	last_id = task1.elements[task1.size - 1];

	printf("=== Setting initial conditions\n");
	init_task1(&device, &task1);

	printf("=== Loading device program\n");
	// Load the device program onto all cores assigned to task 1.
	if (w_load(&device, &task1, "e_main.srec") != E_OK) {
		printf("ERROR: Unable to load device program.\n");
		return 1;
	}

	printf("=== Starting device\n");
	e_start_group(&device);

	printf("=== Starting counter\n");
	mailbox.go = 1;
	w_write(&device, first_id, _MAILBOX_ADDRESS, &mailbox, sizeof(mailbox));

	printf("=== Waiting for last core\n");
	do {
		msleep(100);
		w_read(&device, last_id, _MAILBOX_ADDRESS, &mailbox, sizeof(mailbox));
	} while (mailbox.done == 0);

	printf("Counted to %i (expected %i)\n", mailbox.counter, task1.size - 1);

	printf("=== Finalizing\n");
	w_free_mapper(&mapper);
	w_free_list(&task1);

	e_close(&device);
	e_finalize();

	printf("=== Done\n");
	return 0;
}
