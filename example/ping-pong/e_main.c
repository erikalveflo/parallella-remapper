#include "e_lib.h"

#include "shared.h"

int main(void)
{
	volatile Mailbox *mailbox;
	mailbox = (Mailbox *)_MAILBOX_ADDRESS;

	while (mailbox->ping == 0) {}
	mailbox->pong = 1;

	return 0;
}
