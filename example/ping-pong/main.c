//
// Ping-pong example: the host sends a ping and the device responds with a pong.
//
// In this example a single core is used. This core waits to receive a ping from
// the host, it then responds with a pong. The ping and pong are integers
// located in the cores scratch-pad memory.
//
// The naive assignment algorithm is used in this example. It works well for
// applications that only use a single core.
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

// Allocate a single core for task 1. Use power of two (2, 4, 8...) when adding
// additional tasks. A task is an Epiphany device program.
w_matrix_element_t allocation_matrix[16] = {
	1, 0, 0, 0, // Element 0 corresponds to core (0,0), element 1 to (0,1), and
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

// Initializes and setup all cores assigned to execute task 1.
static void init_task1(e_epiphany_t *device, w_list_t *tasks)
{
	Mailbox mailbox;
	w_core_id_t core_id;
	int i;

	memset(&mailbox, 0, sizeof(mailbox));

	for (i = 0; i < tasks->size; ++i) {
		core_id = tasks->elements[i];

		// Clear the mailbox. Initially, the core's ping and pong value should
		// be zero.
		w_write(device, core_id, _MAILBOX_ADDRESS, &mailbox, sizeof(mailbox));
	}
}

int main(int argc, char *argv[])
{
	e_platform_t platform;
	e_epiphany_t device;

	w_mapper_t mapper;
	w_list_t task1;
	w_core_id_t core_id;

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

	// Use the naive assignment algorithm since we only have a single task.
	// There is also a simulated annealing based approach called w_assign_sa().
	if (w_assign_naive(&mapper) != E_OK) {
		printf("ERROR: Assignment failed.\n");
		return 1;
	}

	w_print_matrix(&mapper.assignment_matrix, "Assignment");
	w_print_matrix(&mapper.mapping_matrix, "Mapping");

	// Find the ID of all cores assigned task 1.
	w_find_in_matrix(&mapper.assignment_matrix, &task1, 1);
	core_id = task1.elements[0];

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

	printf("=== Sending PING\n");
	// Send ping=1 to the core's scratch-pad memory.
	mailbox.ping = 1;
	w_write(&device, core_id, _MAILBOX_ADDRESS, &mailbox, sizeof(mailbox));

	printf("=== Waiting for PONG\n");
	// Read until the core sets pong=1 in its scratch-pad memory.
	do {
		msleep(100);
		w_read(&device, core_id, _MAILBOX_ADDRESS, &mailbox, sizeof(mailbox));
	} while (mailbox.pong == 0);

	printf("PONG received\n");

	printf("=== Finalizing\n");
	w_free_mapper(&mapper);
	w_free_list(&task1);

	e_close(&device);
	e_finalize();

	printf("=== Done\n");
	return 0;
}
