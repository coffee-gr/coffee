# Produced by Ben Whale 2018
#
# The resulting shared library must be copied, manually, to ../../lib

src = $(wildcard src/*.c)
obj = $(src:src/%.c=obj/%.o) 
dep = $(obj:obj/%.o=dep/%.d)

CFLAGS = -fPIC -std=c99 -g -Wall
CINCLUDES = -Iincludes
LDFLAGS = -lhdf5 -shared

-include $(dep)   # include all dep files in the makefile

# rule to generate a dep file by using the C preprocessor
# (see man cpp for details on the -MM and -MT options)
dep/%.d: src/%.c
	@$(CPP) $(CFLAGS) $(CINCLUDES) $< -MM -MT $(@:.d=.o) >$@

libboris.so.1.0.1: $(obj)
	$(CC) -o $@ $^ $(LDFLAGS)

obj/%.o: src/%.c
	$(CC) $(CFLAGS) $(CINCLUDES) -c -o $@ $^

.PHONY: clean
clean:
	rm -f $(obj) myprog

.PHONY: cleandep
cleandep:
	rm -f $(dep)
