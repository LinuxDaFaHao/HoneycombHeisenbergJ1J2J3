#pragma once
#include <stdlib.h>
#include <vector>

size_t GetNumofMps();
void Show(std::vector<size_t> v);
bool ParserBondDimension(int argc, char *argv[],
                         std::vector<size_t>& D_set);
