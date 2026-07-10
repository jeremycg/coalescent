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

// The panel captions live in their own `widget::Widget` subclass (one per module)
// added as a child, rather than drawn inline in ModuleWidget::draw(). Fonts are
// loaded per-frame inside it (see the widget). Call after setPanel() so the
// ModuleWidget's box.size is set.
// (A FramebufferWidget wrapper would cache the unchanging text, but it rendered the
// labels at the wrong scale here, so the labels draw normally each frame.)
template <typename TLabels>
inline void addPanelLabels(ModuleWidget* mw) {
    auto* labels = new TLabels;
    labels->box.size = mw->box.size;
    mw->addChild(labels);
}
