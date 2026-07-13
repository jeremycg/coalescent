#include "plugin.hpp"

Plugin* pluginInstance;

void init(Plugin* p) {
    pluginInstance = p;
    p->addModel(modelGENDYN);
    p->addModel(modelHaptik);
    p->addModel(modelAxon);
    p->addModel(modelSoma);
    p->addModel(modelOperon);
    p->addModel(modelBunnies);
    p->addModel(modelFoxes);
}
