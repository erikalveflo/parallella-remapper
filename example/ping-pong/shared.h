#pragma once

#define _MAILBOX_ADDRESS 0x7000u

typedef struct {
	char ping;
	char pong;
} Mailbox;
