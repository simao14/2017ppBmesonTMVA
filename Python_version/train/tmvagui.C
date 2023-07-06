
#include <string>
#include "TMVA/TMVAGui.h"

int tmvagui(std::string inputname)
{
  TMVA::TMVAGui( inputname.c_str() );
  return 0;
}

