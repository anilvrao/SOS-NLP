ifeq ($(origin FC),default)
FC := gfortran
endif
ifeq ($(origin CC),default)
CC := cc
endif
CMAKE ?= cmake
AR ?= ar
ARFLAGS ?= rcs

FC_NAME := $(notdir $(FC))
ifneq (,$(filter nagfor unagfor,$(FC_NAME)))
FFLAGS ?= -O2 -dusty
else
FFLAGS ?= -O2 -std=legacy -fallow-argument-mismatch
endif

HIGHS_SRC := third_party/HiGHS
HIGHS_BUILD := third_party/highs-build
HIGHS_PREFIX := third_party/highs-install
HIGHS_LIB := $(HIGHS_PREFIX)/lib/libhighs.a

CFLAGS ?= -O2 -I$(HIGHS_PREFIX)/include -I$(HIGHS_PREFIX)/include/highs

SRC_DIR := src/soscode
BRIDGE_DIR := src/highs_bridge
OBJ_DIR ?= build/obj
LIB ?= build/libsos_nlp.a

EXCLUDED_SOURCES := $(SRC_DIR)/lpcore.f $(SRC_DIR)/qpcore.f $(SRC_DIR)/qpopt.f
SOURCES := $(filter-out $(EXCLUDED_SOURCES),$(wildcard $(SRC_DIR)/*.f))
CSOURCES := $(wildcard $(BRIDGE_DIR)/*.c)
OBJECTS := $(patsubst $(SRC_DIR)/%.f,$(OBJ_DIR)/%.o,$(SOURCES))
COBJECTS := $(patsubst $(BRIDGE_DIR)/%.c,$(OBJ_DIR)/%.o,$(CSOURCES))

.PHONY: all highs lib clean distclean manifest

all: lib

highs: $(HIGHS_LIB)

$(HIGHS_LIB):
	$(CMAKE) -S $(HIGHS_SRC) -B $(HIGHS_BUILD) \
		-DCMAKE_BUILD_TYPE=Release \
		-DCMAKE_INSTALL_PREFIX=$(abspath $(HIGHS_PREFIX)) \
		-DBUILD_SHARED_LIBS=OFF \
		-DBUILD_CXX_EXE=OFF \
		-DBUILD_EXAMPLES=OFF \
		-DBUILD_TESTING=OFF \
		-DFAST_BUILD=ON
	$(CMAKE) --build $(HIGHS_BUILD) --target install

lib: $(LIB)

$(LIB): $(HIGHS_LIB) $(OBJECTS) $(COBJECTS) | build
	$(AR) $(ARFLAGS) $@ $(OBJECTS) $(COBJECTS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f | $(OBJ_DIR)
	cd $(SRC_DIR) && $(FC) $(FFLAGS) -c $*.f -o ../../$@

$(OBJ_DIR)/%.o: $(BRIDGE_DIR)/%.c $(HIGHS_LIB) | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

build:
	mkdir -p $@

$(OBJ_DIR):
	mkdir -p $@

manifest:
	find . \
	  -path './.git' -prune -o \
	  -path './build' -prune -o \
	  -path './third_party/highs-build' -prune -o \
	  -path './third_party/highs-install' -prune -o \
	  -type f -print | sort > MANIFEST.generated.txt

clean:
	rm -rf $(OBJ_DIR) $(LIB)

distclean: clean
	rm -rf build $(HIGHS_BUILD) $(HIGHS_PREFIX)
