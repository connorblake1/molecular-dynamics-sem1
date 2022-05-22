//
//  PVector.hpp
//  TrafficSimulator
//
//  Created by Connor Blake on 8/23/21.
//  Copyright © 2021 Connor Blake. All rights reserved.
//

#ifndef PVector_hpp
#define PVector_hpp
#include <math.h>
#include <stdio.h>
#include <iostream>
namespace phys {
class PVector {
private:
    float x,y,z;
public:
    PVector() {}
    PVector(float ix, float iy, float iz) {
        this->x=ix;
        this->y=iy;
        this->z=iz;}
    PVector(float ix, float iy) {
        this->x=ix;
        this->y=iy;
        this->z=0;
    }
    static PVector * rand2DScale(float scale);
    static PVector * rand3DScale(float scale);
    void normalize();
    void add(PVector * other);
    float getX() {return this->x;}
    float getY() {return this->y;}
    float getZ() {return this->z;}
    void setX(float inx) {this->x=inx;}
    void setY(float iny) {this->y=iny;}
    void setZ(float inz) {this->z=inz;}
    void printIt() {std::cout << x << " " << y << " " << z << std::endl;}
    float mag();
    float magsq();
    float dot(PVector * other) {return this->x*(*other).getX()+this->y*(*other).getY()+this->z*(*other).getZ();}
    float theta() {return atan2(this->y,this->x);}
    float phi() {return acos(this->z/sqrt(this->x*this->x+this->y*this->y+this->z*this->z));}
    PVector * cross(PVector * other);
    static PVector * scalm(PVector * other, float scale);
    static PVector * vm(PVector * other, PVector * scaler);
    static float fCross(PVector * other1, PVector * other2);
    static PVector * sub(PVector * other1, PVector * other2);
};}

#endif /* PVector_hpp */
