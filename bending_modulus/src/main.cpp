#include "analyze.h"
#include <unistd.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

int main(int argc, char *argv[])
{
  if (argc == 1) {
    printf("Usage: command [-f input_file]\n");
    return 0;
  }

  int opt;
  std::string input;
  while ((opt = getopt(argc, argv, "f:")) != -1) {
    switch (opt) {
      case 'f':
        input += optarg;
        break;
      default:
        printf("Usage: command [-f input_file]\n");
        return 1;
    }
  }
  if (input.empty()) {
    printf("ERROR: input file is not specified.\n");
    return 1;
  }
  std::ifstream ifs(input, std::ios::in);
  if (!ifs) {
    printf("ERROR: %s cannot open.", input.c_str());
    return 1;
  }

  int beg, end;
  std::string line, str;
  std::string psf, dcd, out;
  std::string lipname;
  std::vector<std::string> heads, tails;
  while (getline(ifs, line)) {
    std::stringstream ss(line);
    ss >> str;
    if (str == "psffile") {
      ss >> psf;
    } else if (str == "dcdfile") {
      ss >> dcd;
    } else if (str == "output") {
      ss >> out;
    } else if (str == "lipid") {
      ss >> lipname;
    } else if (str == "begin") {
      ss >> beg;
    } else if (str == "end") {
      ss >> end;
    } else if (str == "heads") {
      while (!ss.eof()) {
        ss >> str;
        heads.push_back(str);
      }
    } else if (str == "tails") {
      while (!ss.eof()) {
        ss >> str;
        tails.push_back(str);
      }
    }
  }

  LipidBilayer lb(lipname, heads, tails);
  lb.do_analysis(psf.c_str(), dcd.c_str(), out.c_str(), beg, end);
  return 0;
}

