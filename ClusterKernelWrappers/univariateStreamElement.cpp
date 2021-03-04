#include "univariateStreamElement.h"

UnivariateStreamElement::UnivariateStreamElement(const Point &pt){
  coordinates_ = pt;
}

Point UnivariateStreamElement::GetMean() {
  return coordinates_;
}
