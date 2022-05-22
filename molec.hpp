//
//  molec.hpp
//  MD1
//
//  Created by Connor Blake on 8/21/21.
//  Copyright © 2021 Connor Blake. All rights reserved.
//

#ifndef molec_hpp
#define molec_hpp

#include <stdio.h>
#include "PVector.hpp"
#include <SFML/Graphics.hpp>
using namespace phys;
namespace chem {

class molec {
private:
    PVector * r;
    PVector * v;
    PVector * f;
    PVector * centers;
    float rad, m, mi, l, gscale, theta, omega, tq;
    float * charges;
    sf::CircleShape * circles;
    sf::RectangleShape * boxes;
    int molecType, centerNum, barNum; //0 is argon, 1 is nitrogen, 2 is water
public:
    molec() {}
    molec(PVector * inr, PVector * inv, float radin, float min, float fscalin, int type); //nitrogen and water
    float getX();
    float getY();
    float getVX();
    float getVY();
    float getFX();
    float getFY();
    float getRad();
    float getLen();
    float getQ();
    float getM();
    float getMI();
    int getCNum();
    float getTheta();
    float getOmega();
    float getTorque();
    float getCharge(int ind);
    PVector * getR();
    PVector * getV();
    PVector * getF();
    PVector * getC(int ind);
    PVector * getCRotate(int ind);
    void addTheta(float in);
    void addOmega(float in);
    void addTorque(float in);
    void setR(PVector * inr);
    void setTheta(float intheta);
    void setV(PVector * inv);
    void setX(float inx);
    void setY(float iny);
    void setZ(float inz);
    void updateCenters();
    void resetA();
    void disp(sf::RenderWindow * wind);
};
}
#endif /* molec_hpp */
