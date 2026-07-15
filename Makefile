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
# Validate the manifest and run the standalone DSP stability replicas plus the
# RK4 equivalence proof. No Rack SDK required; used by CI and local dev alike.
CHECK_CXX ?= g++ -std=c++17 -O2
# `check` is deliberately SDK-free (the CI guardrail runs it before the Rack SDK is
# downloaded). Anything needing Rack headers (simd/*) goes in `check-simd`, which CI
# runs as a separate post-SDK step.
.PHONY: check check-simd
check:
	jq . plugin.json >/dev/null
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
	$(CHECK_CXX) -funsafe-math-optimizations tools/integrator_equiv.cpp -o /tmp/coalescent_check_equiv && /tmp/coalescent_check_equiv

# Needs Rack headers (simd). Run after the SDK is available: make check-simd
check-simd:
	$(CHECK_CXX) -funsafe-math-optimizations -march=nehalem -DARCH_X64 -DARCH_LIN -I$(RACK_DIR)/include -I$(RACK_DIR)/dep/include tools/simd_equiv.cpp -o /tmp/coalescent_check_simd && /tmp/coalescent_check_simd
