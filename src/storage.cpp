#include "storage.h"
#include <bits/stdc++.h>

void add_to_arr(int x, int j, int beta,std::vector<unsigned> &arr){

    unsigned temp = j*beta;
    unsigned idx = (temp) >> 5;
    unsigned rem=temp&(0x0000001F); 

    unsigned add = rem+beta;

    if(add > 32){
        arr[idx] = arr[idx] & ( ((1 << rem) - 1) << (32-rem) );
        arr[idx] = arr[idx] | (x >> (add-32));
        arr[idx+1] = (x << (64-add)) | ( arr[idx+1] & ( (1<<64-add) - 1)   )  ;
    }
    else{
        unsigned mask_left,mask_right;
        mask_left = ( (1<<rem) - 1) << (32-rem);
        mask_right = (1<<(32-add)) -1 ;
        arr[idx] = (arr[idx] & mask_left) | (arr[idx] & mask_right);
        arr[idx] = arr[idx] | (x<<(32-add)) ;
    } 
}


int extract_from_arr(int j,int beta,std::vector<unsigned> &arr){
    unsigned temp = j*beta;
    unsigned idx = (temp) >> 5;
    unsigned rem=temp&(0x0000001F);
    
    int res=0;
    
    unsigned add = rem + beta;

    if(add > 32){
        unsigned  mask=0;
        mask = (1 << (32-rem)) - 1;
        int lastXbits = arr[idx] & mask;
        res = (lastXbits << (add-32)) | (arr[idx+1] >> (64-add));       
    }
    else{
        res = (arr[idx]<<rem)>>(32-beta);
    }
    return res;
}

std::vector<int> get_orig_arr(int beta,std::vector<unsigned> &arr,int Nelements){

    std::vector<int> v;
    for(int i=0;i<Nelements;i++){
        int temp = extract_from_arr(i,beta,arr);
        v.push_back(temp);
    }
    return v;
}
