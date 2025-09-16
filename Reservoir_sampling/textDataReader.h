//
// Created by tomev on 15/03/2021.
//

#ifndef DEDSTA_TEXTDATAREADER_H
#define DEDSTA_TEXTDATAREADER_H

#include <fstream>
#include "dataReader.h"


using std::string;
using std::vector;

class TextDataReader : public dataReader {

  // This file is for 1D real data, each value in separate line of the file.
  public:

    explicit TextDataReader(const string &path_to_text_file, const int &dimension = 1);
    ~TextDataReader();

    void getNextRawDatum(void *target) override;
    void gatherAttributesData(void *attributes) override;
    bool hasMoreData() override;

    std::vector<std::string>* getAttributesOrder() override;

  protected:
    vector<string> attributesOrder;
    vector<string> lines;
    int i = 0;
    int dimension_ = 1;
};

#endif //DEDSTA_TEXTDATAREADER_H
