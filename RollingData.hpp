//
//  RollingData.hpp
//  MD1
//
//  Created by Connor Blake on 11/28/21.
//  Copyright © 2021 Connor Blake. All rights reserved.
//

#ifndef RollingData_hpp
#define RollingData_hpp
#include <stdio.h>
#include <iostream>
#include <math.h>
class RollingData {
private:
    float * mem;
    int len;
public:
    RollingData() {}
    RollingData(int inLen) {
        this->len = inLen;
        this->mem = new float[len];
        for (int i = 0; i < len; i++) {
            mem[i] = 0.0;
        }
    }
    void push(float in);
    void purge();
    float avg();
    float stdev();
    float variance();
    float * getDat() {return mem;}
};
#endif /* RollingData_hpp */
