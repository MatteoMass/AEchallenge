#pragma once

#include <algorithm>
#include <math.h>
#include <vector>
#include <stdint.h>
#include <iostream>


class MyIndex {
  size_t n;                             
  const std::vector<uint64_t> &data;    
  int sizeBlock;
  std::vector<int> j;


public:

  MyIndex(const std::vector<uint64_t> &data) : n(data.size()), data(data) {
    int nBlock = log10(data.back()+1)+1;
    j.resize(nBlock+2,-1);
    j[0] = 0;
    j.back() = data.size()-1;

    for(int x = 0; x < data.size(); ++x){
        int m = log10(data[x]+1);
        if (j[m+1] != -1){
            ;
        }
        else{
            j[m+1] = x;
        }
    }
    
    int curr = 0;
    for(int x = 1; x < j.size(); ++x){
        if(j[j.size()-x] == -1){
            j[j.size()-x] = curr;
        }else{
            curr = j[j.size()-x];
        }
    }
  }
  


  uint64_t nextGEQ(uint64_t key) const { 
    uint64_t search = log10(key+1);
    uint64_t start = j[search+1];
    uint64_t end = j[search+2];
    return *std::lower_bound(data.begin()+start, data.begin()+end, key);

  }


  size_t size_in_bytes() const { 
    return sizeof(this) + sizeof(int) * j.size();
    
  }
};