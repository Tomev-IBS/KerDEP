//
// Created by tomev on 15/03/2021.
//

#include <sstream>

#include "textDataReader.h"
#include "../groupingThread/kMedoidsAlgorithm/numericalAttributeData.h"

#include <qDebug>

vector<string> tokenize(string s, string del = " ")
{
  vector<string> tokens = {};
  int start = 0;
  int end = s.find(del);
  while (end != -1) {
    tokens.push_back(s.substr(start, end - start));
    start = end + del.size();
    end = s.find(del, start);
  }
  tokens.push_back(s.substr(start, end - start));
  return tokens;
}

TextDataReader::TextDataReader(const string &path_to_text_file) {
  auto opened_file = std::ifstream(path_to_text_file);

  string line;
  while(not opened_file.eof()) {
    if(std::getline(opened_file, line)) {
      lines.push_back(line);
    }
  }

  opened_file.close();
}

TextDataReader::~TextDataReader(){ }

void TextDataReader::getNextRawDatum(void *target) {
  vector<double> *targetPtr = static_cast<vector<double> *>(target);
  targetPtr->clear();

  auto values = tokenize(lines[i++], "; ");

  for(auto value : values){
    targetPtr->push_back(std::stod(value));
  }
}

void TextDataReader::gatherAttributesData(void *attributes) {
  std::unordered_map<std::string, attributeData *> *attrs_ptr =
    static_cast<std::unordered_map<std::string, attributeData *> *>(attributes);
  std::vector<string> attrsNames = {"Val0"}; // 1D
  attrsNames = {"Val0", "Val1"}; // 2D

  for(auto attrName : attrsNames){
    attributesOrder.push_back(attrName);
    (*attrs_ptr)[attrName] = new numericalAttributeData(attrName);
  }
}

bool TextDataReader::hasMoreData() {
  return i < lines.size();
}

std::vector<std::string> *TextDataReader::getAttributesOrder() {
  return &attributesOrder;
}


