//
// Created by tomev on 15/03/2021.
//

#include <sstream>

#include "textDataReader.h"
#include "../groupingThread/kMedoidsAlgorithm/numericalAttributeData.h"

#include <qDebug>

TextDataReader::TextDataReader(const string &path_to_text_file) {
  opened_file_ = std::ifstream(path_to_text_file);
}

TextDataReader::~TextDataReader(){
  opened_file_.close();
}

void TextDataReader::getNextRawDatum(void *target) {
  vector<double> *targetPtr = static_cast<vector<double> *>(target);
  targetPtr->clear();

  string line;
  if(std::getline(opened_file_, line)){
    targetPtr->push_back(std::stod(line));
  }

}

void TextDataReader::gatherAttributesData(void *attributes) {
  std::unordered_map<std::string, attributeData *> *attrs_ptr =
    static_cast<std::unordered_map<std::string, attributeData *> *>(attributes);

  string attrName = "Val0"; // 1D data

  attributesOrder.push_back(attrName);

  (*attrs_ptr)[attrName] = new numericalAttributeData(attrName);
}

bool TextDataReader::hasMoreData() {
  return opened_file_.eof();
}

std::vector<std::string> *TextDataReader::getAttributesOrder() {
  return &attributesOrder;
}


