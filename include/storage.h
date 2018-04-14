#ifndef STORAGE_H
#define STORAGE_H

#include <bits/stdc++.h>

void add_to_arr(int x, int j, int beta,std::vector<unsigned> &arr);
int extract_from_arr(int j,int beta,std::vector<unsigned> &arr);
std::vector<int> get_orig_arr(int beta,std::vector<unsigned> &arr,int Nelements);

#endif