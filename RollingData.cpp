//
//  RollingData.cpp
//  MD1
//
//  Created by Connor Blake on 11/28/21.
//  Copyright © 2021 Connor Blake. All rights reserved.
//

#include "RollingData.hpp"
void RollingData::push(float in) {
    for (int n = len-1; n > 0; n--)  {
        mem[n] = mem[n-1];}
    mem[0] = in;
}
void RollingData::purge() {
    for (int i = 0; i < len; i++) {
        mem[i] = 0.0;}}
float RollingData::avg() {
    float avg = 0;
    for (int i = 0; i < len; i++) {
        avg+=mem[i];}
    return avg/len;}
float RollingData::stdev() {
    return sqrt(variance());}

float RollingData::variance() {
    float var = 0;
    float avg = this->avg();
    for (int i = 0; i < len; i++) {
        var+=pow(mem[i]-avg,2);}
    return var/(len-1);
}
