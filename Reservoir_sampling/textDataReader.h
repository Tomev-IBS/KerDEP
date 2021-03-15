//
// Created by tomev on 15/03/2021.
//

#ifndef KERDEP_TEXTDATAREADER_H
#define KERDEP_TEXTDATAREADER_H

#include <fstream>
#include "dataReader.h"


using std::string;
using std::vector;

class TextDataReader : public dataReader {

  // This file is for 1D real data, each value in separate line of the file.
  public:

    TextDataReader(const string &path_to_text_file);
    ~TextDataReader();

    void getNextRawDatum(void *target) override;
    void gatherAttributesData(void *attributes) override;
    bool hasMoreData() override;

    virtual std::vector<std::string>* getAttributesOrder() override;

  protected:
    vector<string> attributesOrder;
    std::ifstream opened_file_;
};

#endif //KERDEP_TEXTDATAREADER_H
