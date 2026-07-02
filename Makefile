FC ?= gfortran
CC ?= cc
AR ?= ar
ARFLAGS ?= rcs
HIGHS_PREFIX ?= third_party/highs-install

FC_NAME := $(notdir $(FC))
ifneq (,$(filter nagfor unagfor,$(FC_NAME)))
FFLAGS ?= -O2 -dusty
else
FFLAGS ?= -O2 -std=legacy -fallow-argument-mismatch
endif
CFLAGS ?= -O2 -I$(HIGHS_PREFIX)/include -I$(HIGHS_PREFIX)/include/highs

SRC_DIR := src/soscode
BRIDGE_DIR := src/highs_bridge
OBJ_DIR ?= build/obj
LIB ?= build/libsos_nlp.a

EXCLUDED_SOURCES := $(SRC_DIR)/qpopt.f $(SRC_DIR)/qpcore.f $(SRC_DIR)/lpcore.f
SOURCES := $(filter-out $(EXCLUDED_SOURCES),$(wildcard $(SRC_DIR)/*.f))
CSOURCES := $(wildcard $(BRIDGE_DIR)/*.c)
OBJECTS := $(patsubst $(SRC_DIR)/%.f,$(OBJ_DIR)/%.o,$(SOURCES))
COBJECTS := $(patsubst $(BRIDGE_DIR)/%.c,$(OBJ_DIR)/%.o,$(CSOURCES))

.PHONY: all lib clean distclean

all: lib

lib: $(LIB)

$(LIB): $(OBJECTS) $(COBJECTS) | build
	$(AR) $(ARFLAGS) $@ $(OBJECTS) $(COBJECTS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f | $(OBJ_DIR)
	cd $(SRC_DIR) && $(FC) $(FFLAGS) -c $*.f -o ../../$@

$(OBJ_DIR)/%.o: $(BRIDGE_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

build:
	mkdir -p $@

$(OBJ_DIR):
	mkdir -p $@

clean:
	rm -rf $(OBJ_DIR) $(LIB)

distclean: clean
	rm -rf build
