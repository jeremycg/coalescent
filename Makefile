RACK_DIR ?= $(HOME)/Rack2-SDK/Rack-SDK

FLAGS +=
CFLAGS +=
CXXFLAGS +=
LDFLAGS +=

# Resolve #include "plugin.hpp" from nested source dirs (e.g. src/neuron/).
FLAGS += -Isrc

# Recursive: picks up src/*.cpp and src/neuron/*.cpp.
SOURCES += $(shell find src -name '*.cpp')

DISTRIBUTABLES += res
DISTRIBUTABLES += $(wildcard LICENSE*)
DISTRIBUTABLES += $(wildcard README*)
DISTRIBUTABLES += $(wildcard patches)
DISTRIBUTABLES += $(wildcard presets)

include $(RACK_DIR)/plugin.mk

# MinGW binutils 2.30 ld segfaults on the DWARF debug info in these objects;
# disable debug info for Windows cross-builds (appended after plugin.mk so
# -g0 lands after compile.mk's -g and wins).
ifdef ARCH_WIN
  FLAGS += -g0
endif
