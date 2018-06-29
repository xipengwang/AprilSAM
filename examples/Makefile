CC = gcc

CFLAGS = -g -std=gnu99 -Wall -Wno-unused-parameter -Wno-unused-function -O3 -I..
LDFLAGS = -lpthread ../lib/libaprilsam.a -lm
#LDFLAGS = -laprilsam -lpthread -lm

TARGETS := aprilsam_demo aprilsam_tutorial

.PHONY: all
all: $(TARGETS)

aprilsam_demo: aprilsam_demo.o
	@echo "   [$@]"
	@$(CC) -o $@ $^ $(LDFLAGS)

aprilsam_tutorial: aprilsam_tutorial.o ../lib/libaprilsam.a
	@echo "   [$@]"
	@$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	@echo "   $@"
	@$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean
clean:
	@rm -rf *.o $(TARGETS)
