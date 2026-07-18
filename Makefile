RACK_DIR ?= $(HOME)/Rack2-SDK/Rack-SDK

FLAGS +=
CFLAGS +=
CXXFLAGS +=
LDFLAGS +=

# Resolve #include "plugin.hpp" from nested source dirs (e.g. src/neuron/).
FLAGS += -Isrc

# Recursive: picks up src/*.cpp and src/neuron/*.cpp. Sorted for a reproducible
# link order across machines.
SOURCES += $(sort $(shell find src -name '*.cpp'))

DISTRIBUTABLES += res
DISTRIBUTABLES += $(wildcard LICENSE*)
DISTRIBUTABLES += $(wildcard README*)
DISTRIBUTABLES += $(wildcard patches)
DISTRIBUTABLES += $(wildcard presets)

# The Rack SDK is only needed to build the plugin. Skip its include for the
# standalone `make check` target so the DSP guardrail can run in CI (and on a
# dev box) without downloading the SDK.
ifneq ($(filter check,$(MAKECMDGOALS)),check)
include $(RACK_DIR)/plugin.mk

# MinGW binutils 2.30 ld segfaults on the DWARF debug info in these objects;
# disable debug info for Windows cross-builds (appended after plugin.mk so
# -g0 lands after compile.mk's -g and wins).
ifdef ARCH_WIN
  FLAGS += -g0
endif
endif

# ── Non-Rack guardrail ──────────────────────────────────────────────────────
# Validate the manifest and run contracts against the shared SDK-free production
# DSP cores, plus the analytic RK4 contract. No Rack SDK required; used by CI and
# local dev alike.
CHECK_CXX ?= g++ -std=c++17 -O2
# `check` is deliberately SDK-free (the CI guardrail runs it before the Rack SDK is
# downloaded). Anything needing Rack headers (simd/*) goes in `check-simd`, which CI
# runs as a separate post-SDK step.
.PHONY: check check-simd check-rack
check:
	jq . plugin.json >/dev/null
	python3 tools/validate_assets.py
	python3 tools/check_patch_reproducibility.py
	python3 tools/check_shared_core_usage.py
	$(CHECK_CXX) tools/completed_path_test.cpp -o /tmp/coalescent_check_path && /tmp/coalescent_check_path
	$(CHECK_CXX) tools/stability/gendyn.cpp -o /tmp/coalescent_check_gendyn && /tmp/coalescent_check_gendyn
	$(CHECK_CXX) tools/stability/axon.cpp   -o /tmp/coalescent_check_axon   && /tmp/coalescent_check_axon
	$(CHECK_CXX) tools/stability/soma.cpp   -o /tmp/coalescent_check_soma   && /tmp/coalescent_check_soma
	$(CHECK_CXX) tools/stability/haptik.cpp -o /tmp/coalescent_check_haptik && /tmp/coalescent_check_haptik
	$(CHECK_CXX) tools/stability/operon.cpp -o /tmp/coalescent_check_operon && /tmp/coalescent_check_operon
	$(CHECK_CXX) tools/stability/bunnies.cpp -o /tmp/coalescent_check_bunnies && /tmp/coalescent_check_bunnies
	$(CHECK_CXX) tools/stability/foxes.cpp  -o /tmp/coalescent_check_foxes  && /tmp/coalescent_check_foxes
	$(CHECK_CXX) -std=c++11 tools/stability/finches.cpp -o /tmp/coalescent_check_finches && /tmp/coalescent_check_finches
	$(CHECK_CXX) -std=c++11 tools/stability/islands.cpp -o /tmp/coalescent_check_islands && /tmp/coalescent_check_islands
	$(CHECK_CXX) -std=c++11 tools/stability/archipelago.cpp -o /tmp/coalescent_check_archipelago && /tmp/coalescent_check_archipelago
	$(CHECK_CXX) -std=c++11 tools/stability/lineages.cpp -o /tmp/coalescent_check_lineages && /tmp/coalescent_check_lineages
	$(CHECK_CXX) -funsafe-math-optimizations tools/integrator_equiv.cpp -o /tmp/coalescent_check_equiv && /tmp/coalescent_check_equiv

# Needs Rack headers (simd). Run after the SDK is available: make check-simd
check-simd:
	$(CHECK_CXX) -funsafe-math-optimizations -march=nehalem -DARCH_X64 -DARCH_LIN -I$(RACK_DIR)/include -I$(RACK_DIR)/dep/include tools/simd_equiv.cpp -o /tmp/coalescent_check_simd && /tmp/coalescent_check_simd

# Linux-only integration harnesses compile all eleven production Rack wrappers
# and use Engine::resetModule() for the same ResetEvent path as context-menu
# Initialize.
RACK_CHECK_CXX ?= $(CXX)
RACK_CHECK_FLAGS ?= -std=c++17 -O2 -funsafe-math-optimizations -Isrc -I$(RACK_DIR)/include -I$(RACK_DIR)/dep/include
RACK_CHECK_LDFLAGS ?= -L$(RACK_DIR) -Wl,-rpath,$(abspath $(RACK_DIR)) -Wl,--allow-shlib-undefined -lRack
RACK_CHECK_LIBRARY_PATH = $(abspath $(RACK_DIR))$(if $(RACK_RUNTIME_LIBRARY_PATH),:$(RACK_RUNTIME_LIBRARY_PATH))

ifeq ($(shell uname -s),Linux)
check-rack:
	$(RACK_CHECK_CXX) $(RACK_CHECK_FLAGS) tools/stability/finches_rack.cpp $(RACK_CHECK_LDFLAGS) -o /tmp/coalescent_check_finches_rack
	$(RACK_CHECK_CXX) $(RACK_CHECK_FLAGS) tools/stability/islands_rack.cpp $(RACK_CHECK_LDFLAGS) -o /tmp/coalescent_check_islands_rack
	$(RACK_CHECK_CXX) $(RACK_CHECK_FLAGS) tools/stability/lineages_rack.cpp $(RACK_CHECK_LDFLAGS) -o /tmp/coalescent_check_lineages_rack
	$(RACK_CHECK_CXX) $(RACK_CHECK_FLAGS) tools/stability/wrappers_rack.cpp $(RACK_CHECK_LDFLAGS) -o /tmp/coalescent_check_wrappers_rack
	env LD_LIBRARY_PATH="$(RACK_CHECK_LIBRARY_PATH)$${LD_LIBRARY_PATH:+:$${LD_LIBRARY_PATH}}" /tmp/coalescent_check_finches_rack
	env LD_LIBRARY_PATH="$(RACK_CHECK_LIBRARY_PATH)$${LD_LIBRARY_PATH:+:$${LD_LIBRARY_PATH}}" /tmp/coalescent_check_islands_rack
	env LD_LIBRARY_PATH="$(RACK_CHECK_LIBRARY_PATH)$${LD_LIBRARY_PATH:+:$${LD_LIBRARY_PATH}}" /tmp/coalescent_check_lineages_rack
	env LD_LIBRARY_PATH="$(RACK_CHECK_LIBRARY_PATH)$${LD_LIBRARY_PATH:+:$${LD_LIBRARY_PATH}}" /tmp/coalescent_check_wrappers_rack
else
check-rack:
	@echo "check-rack is available on Linux only"; exit 1
endif
