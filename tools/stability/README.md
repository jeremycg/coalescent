# Stability Test Architecture

Standalone stability tests and auxiliary render, profiling, and equivalence
tools must execute the same SDK-free production core used by the Rack module.
They must not maintain a second implementation of model constants, derivatives,
integration scheduling, state repair, or event logic.

Tests should contribute independent evidence around that implementation:

- analytical identities, conserved quantities, and known equilibria;
- higher-accuracy or differently formulated numerical references;
- state, finiteness, convergence, performance, and musical-contract properties;
- explicit historical mutants that prove a regression test discriminates.

Historical mutants must be labelled as deliberately incorrect and must never be
used as the passing production path. Rack-only lifecycle, parameter, polyphony,
and serialization behavior belongs in a Rack integration harness.

`make check-rack` supplies that integration layer for all eleven wrappers. The
focused Finches, Islands, and Lineages binaries cover their deeper stateful
contracts; the aggregate wrapper binary keeps construction, finite-output,
Initialize, persistence, event, polyphony, and reset behavior covered for the
remaining eight, as applicable, without copying their DSP mathematics.

`tools/check_shared_core_usage.py` records and enforces the production/test core
mapping for all eleven modules and their auxiliary tools. It requires both
sides to include the core and call representative shared entry points, and rejects
known legacy copy patterns. New DSP modules should begin with an SDK-free core
rather than adding another mapping after the fact.
