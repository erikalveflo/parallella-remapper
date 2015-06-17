#include "e_lib.h"

#include "shared.h"

int main(void)
{
	volatile Mailbox *mailbox;
	volatile Mailbox *next;

	mailbox = (Mailbox *)_MAILBOX_ADDRESS;
	next = (Mailbox *)mailbox->next_mailbox;

	while (mailbox->go == 0) {}

	next->counter = mailbox->counter + 1;
	next->go = 1;

	mailbox->done = 1;

	return 0;
}
