#pragma once
#include <rack.hpp>

using namespace rack;

extern Plugin* pluginInstance;

// Fluctuations series — oscillator collection.
extern Model* modelGENDYN;   // stochastic (Xenakis GENDYN)
extern Model* modelHaptik;   // physical-modelling ring
extern Model* modelAxon;     // FitzHugh–Nagumo neuron
extern Model* modelSoma;     // Hindmarsh–Rose neuron
extern Model* modelOperon;   // repressilator gene circuit (three-phase)
