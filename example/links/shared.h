#pragma once

#define _MAILBOX_ADDRESS 0x7000u

typedef struct {
	char counter;
	char go;
	char done;
	unsigned int next_mailbox; // Don't use a Mailbox pointer here.
} Mailbox;
