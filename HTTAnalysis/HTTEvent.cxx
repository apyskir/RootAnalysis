#include "HTTEvent.h"

//////////////////////////////////////////////
//////////////////////////////////////////////
Wtau::Wtau(){
  clear();
}
Wtau::~Wtau(){;}

void Wtau::clear(){

  pt_ = -999;
}

Wmu::Wmu(){
  clear();
}
Wmu::~Wmu(){;}

void Wmu::clear(){

  pt_ = -999;
}

Welectron::Welectron(){
  clear();
}
Welectron::~Welectron(){;}

void Welectron::clear(){

  pt_ = -999;
}

Wevent::Wevent(){run_ = 0; lumi_ = 0; event_ = 0;}
Wevent::~Wevent(){;}

Wpair::Wpair(){;}
Wpair::~Wpair(){;}

Wtriggers::Wtriggers(){;}
Wtriggers::~Wtriggers(){;}
