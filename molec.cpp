//
//  molec.cpp
//  MD1
//
//  Created by Connor Blake on 8/21/21.
//  Copyright © 2021 Connor Blake. All rights reserved.
//

#include "molec.hpp"
namespace chem {

molec::molec(PVector * inr, PVector * inv, float radin, float min, float scalin, int type) { //nitrogen and water
    r = new PVector(inr->getX(), inr->getY());
    v = new PVector(inv->getX(), inv->getY());
    f = new PVector(0,0,0);
    this->rad=radin;
    this->m=min;
    //this->theta = -6.283/4; //((double)
    this-> theta = floor(2*(rand())/((double)RAND_MAX))* 6.2834/2 -6.283/4;
    this->omega = 0;
    this->tq = 0;
    this->molecType=type;//0 is argon, 1 is nitrogen, 2 is water
    this->gscale = scalin;
    switch (molecType) {
        case 0: {
            centerNum = 1;
            m *= .6;
            l = 0;
            centers = new PVector[centerNum];
            centers[0] = PVector(0,0);
            charges = new float[centerNum];
            charges[0] = 1;
            circles = new sf::CircleShape[centerNum];
            circles[0] = sf::CircleShape();
            circles[0].setPosition(inr->getX(),inr->getY());
            circles[0].setRadius(rad);
            circles[0].setFillColor(sf::Color::White);
            mi = 0;
            barNum = 0;
            break;}
        case 1: {
            centerNum = 2;
            m *= centerNum;

            l = .25;
            centers = new PVector[centerNum];
            centers[0] = PVector(l/2,0);
            centers[1] = PVector(-l/2,0);
            charges = new float[centerNum];
            charges[0] = 0;
            charges[1] = 0;
            circles = new sf::CircleShape[centerNum];
            circles[0] = sf::CircleShape();
            circles[0].setRadius(rad);
            circles[0].setFillColor(sf::Color::White);
            circles[1] = sf::CircleShape();
            circles[1].setRadius(rad);
            circles[1].setFillColor(sf::Color::White);
            mi = m*l*l/4;
            break;}
        case 2: {
            centerNum = 3;
            m *= centerNum;

            l = .3;
            float r105 = 105/2.*3.1416/180;
            centers = new PVector[centerNum];
            centers[0] = PVector(l*(cos(r105)-.068),l*sin(r105));
            centers[1] = PVector(l*(cos(r105)-.068),-l*sin(r105));
            centers[2] = PVector(-.068*l,0);
            charges = new float[centerNum];
            charges[0] = .25;
            charges[1] = .25;
            charges[2] = -.5;
            circles = new sf::CircleShape[centerNum];
            circles[0] = sf::CircleShape();
            circles[0].setRadius(rad);
            circles[0].setFillColor(sf::Color::White);
            circles[1] = sf::CircleShape();
            circles[1].setRadius(rad);
            circles[1].setFillColor(sf::Color::White);
            circles[2] = sf::CircleShape();
            circles[2].setRadius(rad);
            circles[2].setFillColor(sf::Color::White);
            mi = m*l*l*.1065;
//            barNum = 1;
//            boxes = new sf::RectangleShape[barNum];
//            boxes[0] = sf::RectangleShape();
//            boxes[0].setSize(sf::Vector2f(l*gscale,l*gscale/4));
//            boxes[0].setPosition(r->getX(), r->getY());
            break;}
        case 3: {
            //straight 9
            centerNum = 5;
            l = .3;
            centers = new PVector[centerNum];
            charges = new float[centerNum];
            circles = new sf::CircleShape[centerNum];
            int pi = 0;
            for (int i = 0; i < centerNum; i++) {
                centers[i] = PVector((-1+(i%2)*2)*((i+1)/2)*l,0);
                pi += ((i+1)/2)*((i+1)/2);
                charges[i] = -10;
                circles[i] = sf::CircleShape();
                circles[i].setRadius(rad);
                circles[i].setFillColor(sf::Color::White);
            }
            charges[centerNum-1] = 1;
            circles[centerNum-1].setFillColor(sf::Color::Yellow);
            mi = pi*m*l*l;
            
        //bent 6
//            centerNum = 6;
//            centers = new PVector[centerNum];
//            charges = new float[centerNum];
//            circles = new sf::CircleShape[centerNum];
//            l = .35;
//            centers[0] = PVector(l,l);
//            centers[1] = PVector(l,0);
//            centers[2] = PVector(l,-l);
//            centers[3] = PVector(0,0);
//            centers[4] = PVector(-l,0);
//            centers[5] = PVector(-2*l,0);
//            charges[0] = 1;
//            charges[1] = 1;
//            charges[2] = 1;
//            charges[3] = -10;
//            charges[4] = -10;
//            charges[5] = -10;
//            circles[0] = sf::CircleShape();
//            circles[0].setRadius(rad);
//            circles[0].setFillColor(sf::Color::Yellow);
//            circles[1] = sf::CircleShape();
//            circles[1].setRadius(rad);
//            circles[1].setFillColor(sf::Color::Yellow);
//            circles[2] = sf::CircleShape();
//            circles[2].setRadius(rad);
//            circles[2].setFillColor(sf::Color::Yellow);
//            circles[3] = sf::CircleShape();
//            circles[3].setRadius(rad);
//            circles[3].setFillColor(sf::Color::White);
//            circles[4] = sf::CircleShape();
//            circles[4].setRadius(rad);
//            circles[4].setFillColor(sf::Color::White);
//            circles[5] = sf::CircleShape();
//            circles[5].setRadius(rad);
//            circles[5].setFillColor(sf::Color::White);
//            mi = 10*l*l*m;

            
        //bent 7
//            centerNum = 7;
//            centers = new PVector[centerNum];
//            charges = new float[centerNum];
//            circles = new sf::CircleShape[centerNum];
//            l = .35;
//            centers[0] = PVector(l*1.428,l);
//            centers[1] = PVector(l*1.428,0);
//            centers[2] = PVector(l*1.428,-l);
//            centers[3] = PVector(l*.428,0);
//            centers[4] = PVector(l*-.571,0);
//            centers[5] = PVector(l*-1.571,0);
//            centers[6] = PVector(l*-2.571,0);
//            charges[0] = 1;
//            charges[1] = 1;
//            charges[2] = 1;
//            charges[3] = -10;
//            charges[4] = -10;
//            charges[5] = -10;
//            charges[6] = -10;
//            circles[0] = sf::CircleShape();
//            circles[0].setRadius(rad);
//            circles[0].setFillColor(sf::Color::Yellow);
//            circles[1] = sf::CircleShape();
//            circles[1].setRadius(rad);
//            circles[1].setFillColor(sf::Color::Yellow);
//            circles[2] = sf::CircleShape();
//            circles[2].setRadius(rad);
//            circles[2].setFillColor(sf::Color::Yellow);
//            circles[3] = sf::CircleShape();
//            circles[3].setRadius(rad);
//            circles[3].setFillColor(sf::Color::White);
//            circles[4] = sf::CircleShape();
//            circles[4].setRadius(rad);
//            circles[4].setFillColor(sf::Color::White);
//            circles[5] = sf::CircleShape();
//            circles[5].setRadius(rad);
//            circles[5].setFillColor(sf::Color::White);
//            circles[6] = sf::CircleShape();
//            circles[6].setRadius(rad);
//            circles[6].setFillColor(sf::Color::White);
//            mi = 17.7*l*l*m;
//
            
        
//            centers[0] = PVector(0,0);
//            centers[1] = PVector(-l,0);
//            centers[2] = PVector(l,0);
//            centers[3] = PVector(2*l,0);
//            centers[4] = PVector(-2*l,0);
//            centers[5] = PVector(3*l,0);
//            centers[6] = PVector(-3*l,0);
//            charges[1] = -1;
//            charges[2] = -1;
//            charges[3] = -1;
//            charges[4] = -1;
//            charges[5] = -1;
//            charges[6] = 6;
//
//            circles[0] = sf::CircleShape();
//            circles[0].setRadius(rad);
//            circles[0].setFillColor(sf::Color::White);
//            circles[1] = sf::CircleShape();
//            circles[1].setRadius(rad);
//            circles[1].setFillColor(sf::Color::White);
//            circles[2] = sf::CircleShape();
//            circles[2].setRadius(rad);
//            circles[2].setFillColor(sf::Color::White);
//            circles[3] = sf::CircleShape();
//            circles[3].setRadius(rad);
//            circles[3].setFillColor(sf::Color::White);
//            circles[4] = sf::CircleShape();
//            circles[4].setRadius(rad);
//            circles[4].setFillColor(sf::Color::White);
//            circles[5] = sf::CircleShape();
//            circles[5].setRadius(rad);
//            circles[5].setFillColor(sf::Color::White);
//            circles[6] = sf::CircleShape();
//            circles[6].setRadius(rad);
//            circles[6].setFillColor(sf::Color::Yellow);
//            //mi = m*l*l*2; //3
//            //mi = 10*m*l*l; //5
//            mi = 18*m*l*l; //7
            m *= centerNum;

        }
            omega = -200 + ((double) rand())/((double)RAND_MAX)* 400;
    }
//    centerNum=2;
//    centers = new PVector[centerNum];
//    centers[0] = *r;
//    circles = new sf::CircleShape[centerNum];
//    circles[0] = sf::CircleShape();
//    circles[0].setPosition(inr->getX(),inr->getY());
//    circles[0].setRadius(rad);
//    circles[0].setFillColor(sf::Color::White);
    
}

float molec::getX() {return this->r->getX();}
float molec::getY() {return this->r->getY();}
float molec::getVX() {return this->v->getX();}
float molec::getVY() {return this->v->getY();}
float molec::getFX() {return this->f->getX();}
float molec::getFY() {return this->f->getY();}
float molec::getRad() {return this->rad;}
float molec::getLen() {return this->l;}
float molec::getM() {return this->m;}
float molec::getMI() {return this->mi;}
int molec::getCNum() {return this->centerNum;}
float molec::getTheta() {return theta;}
float molec::getOmega() {return omega;}
float molec::getTorque() {return tq;}
float molec::getCharge(int ind) {return charges[ind];}
PVector *molec::getR() {return r;}
PVector *molec::getV() {return v;}
PVector *molec::getF() {return f;}
PVector *molec::getC(int ind) {return &centers[ind];}
PVector *molec::getCRotate(int ind) {
    float l = this->centers[ind].mag();
    float sAngle = atan2(centers[ind].getY(), centers[ind].getX());
    return new PVector(l*cos(sAngle+this->theta),l*sin(sAngle+this->theta));
}
void molec::addTheta(float in) {this->theta+=in;}
void molec::addOmega(float in) {this->omega+=in;}
void molec::addTorque(float in) {this->tq+=in;}
void molec::setR(PVector * inr) {r = inr;}
void molec::setTheta(float intheta) {theta=intheta;}
void molec::setV(PVector * inv) {v = inv;}
void molec::setX(float inx) {r->setX(inx);}
void molec::setY(float iny) {r->setY(iny);}
void molec::setZ(float inz) {r->setZ(inz);}
void molec::resetA() {f->setX(0); f->setY(0); f->setZ(0); tq = 0;}
void molec::updateCenters() {
    for (int i = 0; i < centerNum; i++) {
        float cost = cos(theta);
        float sint = sin(theta);
        float cix = centers[i].getX();
        float ciy = centers[i].getY();
        circles[i].setPosition(gscale*(r->getX()+cix*cost-ciy*sint),gscale*( r->getY()+cix*sint+ciy*cost));}
//    switch(molecType) {
//        case 0: {break;}
//        case 1: {
//            boxes[0].setPosition(circles[1].getPosition().x+rad-l*gscale/8, circles[1].getPosition().y+rad-l*gscale/8);
//            boxes[0].setRotation(theta);
//            break;}
//        case 2: {
//            break;}}
    
}
void molec::disp(sf::RenderWindow * wind) {
    updateCenters();
    //for (int i = 0; i < barNum; i++) {wind->draw(boxes[i]);}
    for (int i = 0; i < centerNum; i++) {wind->draw(circles[i]);}}
}
