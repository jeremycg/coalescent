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
extern Model* modelBunnies;  // predator-prey (Lotka-Volterra / Rosenzweig-MacArthur)

// The static panel captions never change, so render them once into a
// FramebufferWidget instead of laying out nvgText every frame. Each module
// defines its labels in a `widget::Widget` subclass; this wires one up, sized to
// the panel. Call after setPanel() (so the ModuleWidget's box.size is set).
template <typename TLabels>
inline void addFramebufferedLabels(ModuleWidget* mw) {
    auto* fb = new widget::FramebufferWidget;
    fb->oversample = 2.f;                 // crisp small text when Rack zooms in
    fb->box.size = mw->box.size;
    auto* labels = new TLabels;
    labels->box.size = mw->box.size;
    fb->addChild(labels);
    mw->addChild(fb);
}
