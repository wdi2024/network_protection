CCOMP = gcc
#CFLAGS = -O4 -DNDEBUG -DEXCESS_TYPE_LONG -DPRINT_STAT -DCHECK_SOLUTION -Wall -lm
CFLAGS = -g -DPRINT_FLOW -DEXCESS_TYPE_LONG -DPRINT_STAT -DCHECK_SOLUTION -Wall -lm

all: proc_trav 
proc_trav: main.c proc_trav.c input_graph.c timer.c
	$(CCOMP) $(CFLAGS) -o main main.c libm.so 
clean: 
	rm -f main *~
