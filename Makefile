PREFIX ?= /usr/local

CC = gcc
AR = ar

CFLAGS = -std=gnu99 -fPIC -Wall -Wno-unused-parameter -Wno-unused-function
CFLAGS += -I. -O3 -fno-strict-overflow -g
LDFLAGS += -lm


APRILSAM_SRCS := $(shell ls ./aprilsam/*.c ./aprilsam/common/*.c ./aprilsam/csparse/*.c)
APRILSAM_HEADERS := $(shell ls ./aprilsam/*.h ./aprilsam/common/*.h ./aprilsam/csparse/*.h)
APRILSAM_OBJS := $(APRILSAM_SRCS:%.c=%.o)

APRILSAM_LIB := ./lib
TARGETS := $(APRILSAM_LIB)/libaprilsam.a $(APRILSAM_LIB)/libaprilsam.so

.PHONY: all
all: $(TARGETS)
	@$(MAKE) -C examples all


.PHONY: install
install: $(APRILSAM_LIB)/libaprilsam.so
		@chmod +x install.sh
		@cp  $(APRILSAM_LIB)/libaprilsam.so $(PREFIX)/lib/
		@./install.sh $(PREFIX)/include $(APRILSAM_HEADERS)
		@sed 's:^prefix=$$:prefix=$(PREFIX):' < aprilsam.pc.in > aprilsam.pc
		@./install.sh $(PREFIX)/lib/pkgconfig aprilsam.pc
		@rm aprilsam.pc
		@ldconfig

$(APRILSAM_LIB)/libaprilsam.a: $(APRILSAM_OBJS)
		@echo "   [$@]"
		@$(AR) -cr $@ $(APRILSAM_OBJS)

$(APRILSAM_LIB)/libaprilsam.so: $(APRILSAM_OBJS)
		@echo "   [$@]"
		@$(CC) -fPIC -shared -o $@ $^

%.o: %.c
		@echo "   $@"
		@$(CC) -o $@ -c $< $(CFLAGS)


.PHONY: clean
clean:
		@rm -rf *.o aprilsam/*.o aprilsam/common/*.o aprilsam/csparse/*.o unity/*.o $(TARGETS)
		@$(MAKE) -C examples clean
