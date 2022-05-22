//
//  PVector.cpp
//  TrafficSimulator
//
//  Created by Connor Blake on 8/23/21.
//  Copyright © 2021 Connor Blake. All rights reserved.
//

#include "PVector.hpp"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
namespace phys {
PVector * PVector::rand2DScale(float scale) {
    float t = ((double) rand())/((double) RAND_MAX)* 6.283;
    float x1 = cos(t)*scale;
    float y1 = sin(t)*scale;
    return new PVector(x1,y1);}

PVector * PVector::rand3DScale(float scale) {
    float t = ((double) rand())/((double) RAND_MAX)* 6.283;
    float p = ((double) rand()/RAND_MAX) * 3.1415;
    float x1 = cos(t)*sin(p)*scale;
    float y1 = sin(t)*sin(p)*scale;
    float z1 = cos(p)*scale;
    return new PVector(x1,y1,z1);}

void PVector::normalize() {
    float r = this->mag();
    this->x /= r;
    this->y /= r;
    this->z /= r;}

PVector * PVector::scalm(PVector * other, float scale) {
    PVector * r = new PVector((*other).getX()*scale,(*other).getY()*scale,(*other).getZ()*scale);
    return r;}

PVector * PVector::vm(PVector * other, PVector * scaler) {
    return new PVector(other->getX()*scaler->getX(),other->getY()*scaler->getY(),other->getZ()*scaler->getZ());}

void PVector::add(PVector * other) {
    //return new PVector(this->x + other->getX(),this->y + other->getY(), this->z + other->getZ());
    this->x += (*other).getX();
    this->y += (*other).getY();
    this->z += (*other).getZ();}

float PVector::mag() {return sqrt(x*x+y*y+z*z);}
float PVector::magsq() {return x*x+y*y+z*z;}

//PVector * PVector::cross(PVector * other) {
//
//}
float PVector::fCross(PVector * other1, PVector * other2) {return other1->getX()*other2->getY()-other1->getY()*other2->getX();}



PVector * PVector::sub(PVector * other1,PVector * other2) {
    return new PVector((*other1).getX()-(*other2).getX(),(*other1).getY()-(*other2).getY(),(*other1).getZ()-(*other2).getZ());}
}

