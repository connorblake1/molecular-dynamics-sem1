#include <SFML/Audio.hpp>
#include <SFML/Graphics.hpp>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "ResourcePath.hpp"
#include "molec.hpp"
#include "PVector.hpp"
#include "RollingData.hpp"
#include "SliderSFML.hpp"
using namespace phys;
using namespace chem;
unsigned long t;
constexpr float ts = .00005;
const int WINX = 1000;
const int WINY = 1000;
constexpr float mass = 1;
float scrnl = 50;
float effX = 1;
float effSize = -1;
constexpr float epsilon = 500;
constexpr float sigma = 25;
constexpr float fixd = 1.224;//*sigma;
float k = 900.; //TODO match up the units
float DENSITY = .25;
constexpr float DENSITY2 = 2;
constexpr float rl_argon = .00000000034; //in meters
constexpr float re_argon = .000000000000000000001657; //in units of J/atom
constexpr float rm_argon = .0000000000000000000000000669; // in kg
constexpr float rt_argon = .000000000002161;
constexpr float ert_argon = rt_argon*ts;
float ke= 0;
float pe=  0;
float virsum = 0;
//<yoshida constants>
const float c[4] = {.6756,-.1756,-.1756,.6756};
const float d[3] = {1.3512,-1.7024,1.3512};
//</yoshida constants>
int STATE = -1; //-1 is menu, 0 is simulation
int INTEGRATOR = 0; //0 is leapfrog, 1 is yoshida
int POTENTIAL = 0; //0 is no attraction, 1 is LJ attraction
float TEMPERATURE = 125; // in kelvin
float TEMPERATURE2 = 1;
int NUM = 200;//floor(); //NUMber of molecules
int MOLECULE = 0; //0 is argon, 1 is N2, 2 is water, 3 is phospholipid
int SHOW_VARS = 0; //0 is no vars, 1 is vars (pressure, temperature, energy loss etc)
int SIMULATION = 0; //0 is gas, 1 is phase boundary, 2 is phase transition
int BOUNDARY = 0; //0 is wraparound, 1 is hard box
int TRACER = 0;
bool primed = false; // false if any one of the conditions above is -1, true if ready to go
int avgCycle = 1000;
RollingData eAvg(avgCycle);
RollingData pAvg(avgCycle);
RollingData tAvg(avgCycle);
int tracerCycle = 200;
int tracerIndex = -1;
RollingData tracerX(tracerCycle);
RollingData tracerY(tracerCycle);
bool startAvg = false;
std::stringstream estream, mstream, pstream, tstream;

void initializeMolecules(molec * m) {
    switch (SIMULATION) {
        case -1:{
            break;}
        case 0:{ //liquid
            float scaleV = sqrt(2.*(1.0-1./NUM)*TEMPERATURE);
            std::cout << "scalv " << scaleV << std::endl;
            //effX = 40; //pow(NUM,1./3.)*1.218+1;//liquid argon only
            int unitCell = ceil(sqrt(NUM));
            effX = unitCell/sqrt(DENSITY);
            scrnl = WINX/effX;
            float gap = effX/(unitCell*1.);
            //std::cout << scaleV << " " << effX << " " << scrnl << " " << (effX-1)/perEdge << std::endl;
            effSize = fixd/4*scrnl;
            effSize = fixd/8*scrnl;
            int n = 0;
            std::cout << NUM << std::endl;
            PVector * vs = new PVector[NUM];
            for (int i = 0; i < NUM/2; i++) {
                float t = ((double) rand())/((double)RAND_MAX)* 6.283;
                float x1 = cos(t)*scaleV;
                float y1 = sin(t)*scaleV;
                vs[i] = PVector(x1,y1);
                vs[NUM-1-i] = PVector(-x1,-y1);}
            std::cout << " unit cell: " << unitCell << std::endl;

            for (int xi = 0; xi < unitCell; xi++) {
                for (int yi = 0; yi < unitCell; yi++) {
                    if (n < NUM) {
                        PVector * rin;
                        PVector * vin;
                        if (MOLECULE == 3) {
                            float x = gap*(xi+.5);
                            float y = gap*(yi+.5);
                            int mole = 0;
                            if (n % 10 == 0) {
                                mole = 3;
                                int f = 2;
                                yi+= f;
                                y+=gap*f*.5;}
                            rin = new PVector(x,y);
                            vin = &vs[n];
                            m[n] = molec(rin, vin, effSize, mass,scrnl,mole);}

//                        int x1 = 7;
//                        int i1 = 3;
//                        if (MOLECULE == 3) {
//                            float x = gap*(xi+.5);
//                            float y = gap*(yi+.5);
//                            int mole = 0;
//                            if (yi % x1 == 0 && ((xi+1) % i1 == 0 || (xi-1) % i1 == 0)) {continue;}
//                            else if (xi % i1 == 0 && yi % x1 == 0) {
//                                mole = 3;
//                                int f = 2;
//                                yi+= f;
//                                y+=gap*f*.35;}
//                            rin = new PVector(x,y);
//                            vin = &vs[n];
//                            m[n] = molec(rin, vin, effSize, mass,scrnl,mole);}
                        else {
                            float x = gap*(xi+.5);
                            float y = gap*(yi+.5);
                            rin = new PVector(x,y);
                            vin = &vs[n];
                            m[n] = molec(rin, vin, effSize,mass,scrnl,MOLECULE);
                            }
                        n++;
                        delete rin;
                        
                    }}}
            if (n < NUM-1) {
                NUM = n;
            }
            delete [] vs;
            break;}
        case 1: { //phase boundary
//            float rat = pow(DENSITY/DENSITY2,2);
//            int NUM1 = rat*NUM; // num 1 is the hot gas
//            int NUM2 = NUM-NUM1; //cold majority
//            float scaleV1 = sqrt(2.*(1.0-1./NUM1)*TEMPERATURE); //hot
//            float scaleV2 = sqrt(2.*(1.0-1./NUM2)*TEMPERATURE2); //cold
//
//            int unitCell = ceil(sqrt(2*NUM));
//            effX = unitCell/sqrt(DENSITY);
//            scrnl = WINX/effX;
//            float gap = effX/(unitCell*1.);
//            effSize = 20; //fixd/2*scrnl;
//            std::cout << NUM << " " << NUM1 << " " << NUM2 << std::endl;
//            int n = 0;
//            for (float xi = 0; xi < unitCell; xi++) {
//                for (float yi = unitCell/2; yi < unitCell; yi++) {
//                    if (n < NUM2) {
//                        float x = gap*(xi+.5);
//                        float y = gap*(yi+.5);
//                        PVector * rin = new PVector(x,y);
//                        PVector * vin = PVector::rand2DScale(scaleV2);
//                        PVector * ain = PVector::rand2DScale(0);
//                        m[n] = molec(rin, vin, ain,effSize,1,mass);
//                        n++;
//                        delete rin;
//                        delete vin;
//                        delete ain;}}}
//            for (float xi = 0; xi < unitCell; xi++) {
//                for (float yi = 0; yi < unitCell/2; yi++) {
//                    if (n < NUM) {
//                        float x = gap*(xi+.5);
//                        float y = gap*(yi+.5);
//                        PVector * rin = new PVector(x,y);
//                        PVector * vin = PVector::rand2DScale(scaleV1);
//                        PVector * ain = PVector::rand2DScale(0);
//                        m[n] = molec(rin, vin, ain,effSize,1,mass);
//                        n++;
//                        delete rin;
//                        delete vin;
//                        delete ain;}}}
            break;}
        case 2: {//phase transition
            break;}
    }
}

void showTracer(sf::RenderWindow * wind) {
    for (int i = 1; i < tracerCycle; i++) {
        if (tracerX.getDat()[i] == 0.0) {return;}
        sf::Vertex line[] = { sf::Vertex(sf::Vector2f(tracerX.getDat()[i-1],tracerY.getDat()[i-1]),sf::Color::Red), sf::Vertex(sf::Vector2f(tracerX.getDat()[i],tracerY.getDat()[i]),sf::Color::Red)};
        wind->draw(line, 2, sf::Lines);}}

void ljpotforce(molec * other1, molec * other2, float factor, int i1, int i2) {
    PVector * i1v = other1->getCRotate(i1);
    PVector * i2v = other2->getCRotate(i2);
    float x1 = (*other1).getX()+i1v->getX();
    float x2 = (*other2).getX()+i2v->getX();
    float y1 = (*other1).getY()+i1v->getY();
    float y2 = (*other2).getY()+i2v->getY();
    delete i1v;
    delete i2v;
    float dx = abs(x1-x2);
    float dy = abs(y1-y2);
    if (dx > effX/2) {dx-=effX;}
    if (dy > effX/2) {dy-=effX;}
    if (dx >= fixd*factor || dy >= fixd*factor) {
        return;}
    
    float distX = effX;
    float distY = effX;
    int ix = -1;
    int iy = -1;
    for (int i = -1; i <= 1; i++) {
        float nDistY = abs(y1-y2+effX*i);
        float nDistX = abs(x1-x2+effX*i);
        if (nDistX < distX) {
            distX=nDistX;
            ix=i;}
        if (nDistY < distY) {
            distY=nDistY;
            iy=i;}}
    float n  = sqrt(distX*distX+distY*distY);
    if (n < factor*fixd) {
        float angle = atan2(y1-y2+iy*effX,x1-x2+ix*effX);
        PVector * r = new PVector(n*cos(angle),n*sin(angle));
        
        float fn = 48*(pow(n,-14)-.5*pow(n,-8));
        fn = fn/abs(fn)*fmin(abs(fn), 30000000);
        
//        float qn = 0;
//        if (other1->getCharge(i1) > -5 && other2->getCharge(i2) > -5) {
//            qn = -k*other1->getCharge(i1)*other2->getCharge(i2)/n/n/n;
//            if (qn != 0) {
//                qn = qn/abs(qn)*fmin(abs(qn), 10000000);}}
        
        
        float qn = -k*other1->getCharge(i1)*other2->getCharge(i2)/n/n/n;
        if (qn > 0) {
            qn = qn/abs(qn)*fmin(abs(qn), 10000000);}
        else {
            qn = 0;}
        
        
        PVector * r2 = PVector::scalm(r,fn);
        (*(*other1).getF()).add(r2);
        //std::cout << qn << " " << fn << std::endl;
        PVector * r3 = PVector::scalm(r,qn);
        (*(*other1).getF()).add(r3);
        
        float t1 = PVector::fCross(other1->getCRotate(i1), r2);
        other1->addTorque(t1);
        
        PVector * r1 = PVector::scalm(r2,-1);
        (*(*other2).getF()).add(r1);
        PVector * r4 = PVector::scalm(r3,-1);
        (*(*other2).getF()).add(r4);


        float t2 = PVector::fCross(other2->getCRotate(i2), r1);
        other2->addTorque(t2);
        delete r4;
        delete r3;
        delete r2;
        delete r1;
        delete r;
        
        virsum += fn*n;
        pe += 4*(pow(n,-12)-pow(n,-6));
        if (factor == 1.) {pe += 1;}
    }
}


void generalizedForce(molec * other1, molec * other2) {
   // std::cout << "Cycle" << std::endl;;
    for (int i1 = 0; i1 < other1->getCNum(); i1++) {
        for (int i2 = 0; i2 < other2->getCNum(); i2++) {
            ljpotforce(other1, other2,3*POTENTIAL+1, i1,i2);}}}

float actualDistance(float x1, float y1, float x2, float y2) {
    float distX = effX;
    float distY = effX;
    int ix = -1;
    int iy = -1;
    for (int i = -1; i <= 1; i++) {
        float nDistY = abs(y1-y2+effX*i);
        float nDistX = abs(x1-x2+effX*i);
        if (nDistX < distX) {
            distX=nDistX;}
        if (nDistY < distY) {
            distY=nDistY;}}
    return sqrt(distX*distX+distY*distY);
}


void generalizedCalculation(molec * m) {
    switch (SHOW_VARS) {
        case -1: {
            break;}
        case 0: {
            estream.str(std::string());
            estream.clear();
            mstream.str(std::string());
            mstream.clear();
            tstream.str(std::string());
            tstream.clear();
            pstream.str(std::string());
            pstream.clear();
            estream << "Energy (per molec): " << std::fixed << std::setprecision(2) << (ke+pe)/NUM;
            float mx = 0;
            float my = 0;
            for (int i = 0; i < NUM; i++) {
                mx += m[i].getV()->getX();
                my += m[i].getV()->getY();}
            
            mstream << "Momentum (x,y): (" << std::fixed << std::setprecision(2) << mx << ", " << std::fixed << std::setprecision(2) << my << ")";
            
            tstream << "Temperature: " << std::fixed << std::setprecision(2) << ke/NUM;
            float p = DENSITY*(ke*2+virsum)/NUM/2;//(virsum/2+NUM*TEMPERATURE)/WINX/WINY;
            pstream << "Pressure: " << std::fixed << std::setprecision(2) << p;
            pAvg.push(p);
            eAvg.push(((ke+pe)/NUM));
            tAvg.push(ke/NUM);
            if (startAvg) {
                estream << "\t Avg: " << std::fixed << std::setprecision(3) << eAvg.avg();
                pstream << "\t Avg: " << std::fixed << std::setprecision(3) << pAvg.avg();
                tstream << "\t Avg: " << std::fixed << std::setprecision(3) << tAvg.avg();}
            else {
                if (t > avgCycle) {startAvg = true;}}
            
            break;}
        case 1: {
            break;}}}

void generalizedBoundary(molec * m, int index) {
    bool crosser = false;
    switch (BOUNDARY) {
        case -1: {
            break;}
        case 0: { //tiled
            if (m->getX() > effX) {
                crosser = true;
                m->setX(m->getX() - effX);}
            else if (m->getX() < 0) {
                crosser = true;
                m->setX(m->getX() + effX);}
            if (m->getY() > effX) {
                crosser = true;
                m->setY(m->getY() - effX);}
            else if (m->getY() < 0) {
                crosser = true;
                m->setY(m->getY() + effX);}
            break;}
        case 1: { //bouncy box
            if (m->getX() > effX || m->getX() < 0) {
                m->getV()->setX(m->getV()->getX()*-1.);
            }
            if (m->getY() > effX || m->getY() < 0) {
                m->getV()->setY(m->getV()->getY()*-1.);
            }
            break;}
        case 2: { //layer box
            if (m->getX() > effX) {
                crosser = true;
                m->setX(m->getX() - effX);}
            else if (m->getX() < 0) {
                crosser = true;
                m->setX(m->getX() + effX);}
            if (m->getY() > effX || m->getY() < 0) {
                m->getV()->setY(m->getV()->getY()*-1.);}
            break;}}
        if (TRACER == 0) {
            if (index == tracerIndex && crosser) {
                tracerX.purge();
                tracerY.purge();}}
    }
void generalizedIntegrator(molec * m) {
    switch (INTEGRATOR){
        case -1:
            break;
        case 0: { //leapfrog
            for (int i = 0; i < NUM; i++) {
                for (int j = i+1; j < NUM; j++) {
                    generalizedForce(&m[i],&m[j]);}}
            for (int i = 0; i < NUM; i++) {
                PVector * a1 = PVector::scalm(m[i].getF(),ts/2/m[i].getM());
                (*m[i].getV()).add(a1);
                delete a1;
                PVector * a2 = PVector::scalm(m[i].getV(),(ts));
                (*m[i].getR()).add(a2);
                delete a2;
                if (m[i].getCNum() > 1) {
                    m[i].addOmega(m[i].getTorque()*ts/2/m[i].getMI());
                    m[i].addTheta(m[i].getOmega()*ts);}
                m[i].resetA();}
            virsum = 0;
            for (int i = 0; i < NUM; i++) {
                for (int j = i+1; j < NUM; j++) {
                    generalizedForce(&m[i],&m[j]);}}
            for (int i = 0; i < NUM; i++) {
                PVector * a3 = PVector::scalm(m[i].getF(),ts/2/m[i].getM());
                (*m[i].getV()).add(a3);
                delete a3;
                if (m[i].getCNum() > 1) {
                    m[i].addOmega(m[i].getTorque()*ts/2/m[i].getMI());}
                m[i].resetA();
                generalizedBoundary(&m[i],i);
            }
            break;}
        case 1: { //yoshida
            for (int i = 0; i < NUM; i++) {
                PVector * a1 = PVector::scalm(m[i].getV(),ts*c[0]);
                (*m[i].getR()).add(a1);
                delete a1;
                if (m[i].getCNum() > 1) {m[i].addTheta(m[i].getOmega()*ts*c[0]);}}
            for (int l = 0; l <= 2; l++) {
                if (l==2) {virsum = 0;}
                for (int i = 0; i < NUM; i++) {
                    for (int j = i+1; j < NUM; j++) {
                        generalizedForce(&m[i],&m[j]);}}
                for (int i = 0; i < NUM; i++) {
                    PVector * a2 = PVector::scalm(m[i].getF(),ts*d[l]/m[i].getM());
                    (*m[i].getV()).add(a2);
                    delete a2;
                    PVector * a3 = PVector::scalm(m[i].getV(),ts*c[l+1]);
                    (*m[i].getR()).add(a3);
                    delete a3;
                    if (m[i].getCNum() > 1) {
                        m[i].addOmega(m[i].getTorque()*ts*d[l]/m[i].getMI());
                        m[i].addTheta(m[i].getOmega()*ts*c[l+1]);}
                    m[i].resetA();
                    if (l == 2) {generalizedBoundary(&m[i],i);}}}
            break;}}}



int main(int, char const**)
{
    std::cout << "This is the IS branch." << std::endl;
    //<text initialization>
    int fontsize = 32;
    float lalign = 50;
    float ualign = 100;
    float vspacer = 65;
    float optionOffset = 300;
    sf::Color dullWhite(255,255,255,128);
    sf::Font fontMain;
    fontMain.loadFromFile("/Users/connorblake/Documents/SORTOld/RampCurves/RampCurves/arial.ttf");
    sf::Text title("Molecular Dynamics Simulator 1.1",fontMain,1.5*fontsize);
    title.setFillColor(sf::Color::Blue);
    title.setPosition(lalign, ualign);
    
    sf::Text options("Options", fontMain, fontsize*1.5);
    options.setFillColor(sf::Color::Blue);
    options.setPosition(lalign, ualign+vspacer);
    //integrator
    sf::Text integrator("Integrator",fontMain,fontsize);
    integrator.setFillColor(sf::Color::White);
    integrator.setPosition(lalign,ualign+2*vspacer);
    sf::Text integrator1("Leapfrog (speed)",fontMain,fontsize);
    integrator1.setFillColor(dullWhite);
    integrator1.setPosition(lalign+1*optionOffset,ualign+2*vspacer);
    sf::Text integrator2("Yoshida (accuracy)",fontMain,fontsize);
    integrator2.setFillColor(dullWhite);
    integrator2.setPosition(lalign+2*optionOffset,ualign+2*vspacer);
    //boundary
    sf::Text boundary("Boundary",fontMain,fontsize);
    boundary.setFillColor(sf::Color::White);
    boundary.setPosition(lalign,ualign+3*vspacer);
    sf::Text boundary1("Wraparound",fontMain,fontsize);
    boundary1.setFillColor(dullWhite);
    boundary1.setPosition(lalign+1*optionOffset,ualign+3*vspacer);
    sf::Text boundary2("Confined",fontMain,fontsize);
    boundary2.setFillColor(dullWhite);
    boundary2.setPosition(lalign+5./3.*optionOffset,ualign+3*vspacer);
    sf::Text boundary3("Layer",fontMain,fontsize);
    boundary3.setFillColor(dullWhite);
    boundary3.setPosition(lalign+7./3.*optionOffset,ualign+3*vspacer);
    //Potential
    sf::Text potential("Potential",fontMain,fontsize);
    potential.setFillColor(sf::Color::White);
    potential.setPosition(lalign,ualign+4*vspacer);
    sf::Text potential1("Lennard-Jones (no attraction)",fontMain,fontsize/1.5);
    potential1.setFillColor(dullWhite);
    potential1.setPosition(lalign+1*optionOffset,ualign+4*vspacer);
    sf::Text potential2("Lennard-Jones (attraction)",fontMain,fontsize/1.5);
    potential2.setFillColor(dullWhite);
    potential2.setPosition(lalign+2*optionOffset,ualign+4*vspacer);
    //temperature
    sf::Text temperature("Temperature",fontMain,fontsize);
    temperature.setFillColor(sf::Color::White);
    temperature.setPosition(lalign,ualign+5*vspacer);
    SliderSFML sliderT(lalign+optionOffset,ualign+5*vspacer);
    sliderT.create(1, 500,true);
    //number
    sf::Text number("Number",fontMain,fontsize);
    number.setFillColor(sf::Color::White);
    number.setPosition(lalign,ualign+6*vspacer);
    SliderSFML sliderN(lalign+optionOffset,ualign+6*vspacer);
    sliderN.create(10, 400, true);
    //density
    sf::Text density("Density",fontMain,fontsize);
    density.setFillColor(sf::Color::White);
    density.setPosition(lalign,ualign+7*vspacer);
    SliderSFML sliderD(lalign+optionOffset,ualign+7*vspacer);
    sliderD.create(.25, 5.0,false);
    //Potential
    sf::Text display("Display",fontMain,fontsize);
    display.setFillColor(sf::Color::White);
    display.setPosition(lalign,ualign+8*vspacer);
    sf::Text display1("All Variables",fontMain,fontsize);
    display1.setFillColor(dullWhite);
    display1.setPosition(lalign+1*optionOffset,ualign+8*vspacer);
    sf::Text display2("No Variables",fontMain,fontsize);
    display2.setFillColor(dullWhite);
    display2.setPosition(lalign+2*optionOffset,ualign+8*vspacer);
    //Simulation
    sf::Text simulation("Simulation",fontMain,fontsize);
    simulation.setFillColor(sf::Color::White);
    simulation.setPosition(lalign,ualign+9*vspacer);
    sf::Text simulation1("Liquid",fontMain,fontsize/1.5);
    simulation1.setFillColor(dullWhite);
    simulation1.setPosition(lalign+1*optionOffset,ualign+9*vspacer);
    sf::Text simulation2("Phase Boundary",fontMain,fontsize/1.5);
    simulation2.setFillColor(dullWhite);
    simulation2.setPosition(lalign+5./3.*optionOffset,ualign+9*vspacer);
    sf::Text simulation3("Phase Transition",fontMain,fontsize/1.5);
    simulation3.setFillColor(dullWhite);
    simulation3.setPosition(lalign+7./3.*optionOffset,ualign+9*vspacer);
    //molecule
    sf::Text molecule("Molecule",fontMain,fontsize);
    molecule.setFillColor(sf::Color::White);
    molecule.setPosition(lalign,ualign+10*vspacer);
    sf::Text molecule1("Argon",fontMain,fontsize/2);
    molecule1.setFillColor(dullWhite);
    molecule1.setPosition(lalign+1.*optionOffset,ualign+10*vspacer);
    sf::Text molecule2("Nitrogen",fontMain,fontsize/2);
    molecule2.setFillColor(dullWhite);
    molecule2.setPosition(lalign+4./3.*optionOffset,ualign+10*vspacer);
    sf::Text molecule3("Water",fontMain,fontsize/2);
    molecule3.setFillColor(dullWhite);
    molecule3.setPosition(lalign+6./3.*optionOffset,ualign+10*vspacer);
    sf::Text molecule4("Phospholipid",fontMain,fontsize/2);
    molecule4.setFillColor(dullWhite);
    molecule4.setPosition(lalign+8./3.*optionOffset,ualign+10*vspacer);
    //tracer
    sf::Text tracer("Tracer",fontMain,fontsize);
    tracer.setFillColor(sf::Color::White);
    tracer.setPosition(lalign,ualign+11*vspacer);
    sf::Text tracer1("On",fontMain,fontsize);
    tracer1.setFillColor(dullWhite);
    tracer1.setPosition(lalign+1*optionOffset,ualign+11*vspacer);
    sf::Text tracer2("Off",fontMain,fontsize);
    tracer2.setFillColor(dullWhite);
    tracer2.setPosition(lalign+2*optionOffset,ualign+11*vspacer);
    
    //showvars
    int leftalign = 50;
    int upperalign = 50;
    sf::Text energy("Energy: ", fontMain,fontsize/2);
    energy.setFillColor(sf::Color::White);
    energy.setPosition(leftalign,upperalign);
    sf::Text momentum("Momentum: ", fontMain,fontsize/2);
    momentum.setFillColor(sf::Color::White);
    momentum.setPosition(leftalign,upperalign+.5*vspacer);
    sf::Text temp("Temperature: ", fontMain,fontsize/2);
    temp.setFillColor(sf::Color::White);
    temp.setPosition(leftalign,upperalign+1*vspacer);
    sf::Text pressure("Pressure: ", fontMain,fontsize/2);
    pressure.setFillColor(sf::Color::White);
    pressure.setPosition(leftalign,upperalign+1.5*vspacer);
    
    sf::Text initializeAll("Simulate", fontMain, fontsize);
    initializeAll.setFillColor(sf::Color::Blue);
    initializeAll.setPosition(lalign,ualign+12*vspacer);
    //</text initialization>
    // Create the main window
    sf::RenderWindow window(sf::VideoMode(WINX, WINY), "MD Sim 1");
    // Set the Icon
//    sf::Image icon;
//    if (!icon.loadFromFile("/Users/connorblake/Desktop/MD1-feryaccyafhmklbkcnwifmbynyet/Build/Products/Debug/../../../../../Documents/SORTOld/MD1/MD1/MD_Icon.png")) {
//        return EXIT_FAILURE;
//    }
//    window.setIcon(icon.getSize().x, icon.getSize().y, icon.getPixelsPtr());
    molec* gas;// = new molec[NUM];
    srand( (unsigned)time( NULL ) );
    // Start the game loop
    while (window.isOpen())
    {
        // Process events
        sf::Event event;
        while (window.pollEvent(event))
        {
            //close
            if (event.type == sf::Event::Closed) {
                window.close();}
            // Escape pressed: exit
            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape) {
                window.close();}
            switch (STATE) {
                case -1: {
                    if (event.type == sf::Event::MouseButtonReleased) {
                        sf::Mouse mouse;
                        sf::Vector2i mousey = sf::Mouse::getPosition(window);
                        if (integrator1.getGlobalBounds().contains(mousey.x,mousey.y)) {INTEGRATOR = 0;}
                        else if (integrator2.getGlobalBounds().contains(mousey.x,mousey.y)) {INTEGRATOR = 1;}
                        else if (boundary1.getGlobalBounds().contains(mousey.x,mousey.y)) {BOUNDARY = 0;}
                        else if (boundary2.getGlobalBounds().contains(mousey.x,mousey.y)) {BOUNDARY = 1;}
                        else if (boundary3.getGlobalBounds().contains(mousey.x,mousey.y)) {BOUNDARY = 2;}
                        else if (potential1.getGlobalBounds().contains(mousey.x,mousey.y)) {POTENTIAL = 0;}
                        else if (potential2.getGlobalBounds().contains(mousey.x,mousey.y)) {POTENTIAL = 1;}
                        else if (display1.getGlobalBounds().contains(mousey.x,mousey.y)) {SHOW_VARS = 0;}
                        else if (display2.getGlobalBounds().contains(mousey.x,mousey.y)) {SHOW_VARS = 1;}
                        else if (simulation1.getGlobalBounds().contains(mousey.x,mousey.y)) {SIMULATION = 0;}
                        else if (simulation2.getGlobalBounds().contains(mousey.x,mousey.y)) {SIMULATION = 1;}
                        else if (simulation3.getGlobalBounds().contains(mousey.x,mousey.y)) {SIMULATION = 2;}
                        else if (molecule1.getGlobalBounds().contains(mousey.x,mousey.y)) {MOLECULE = 0;}
                        else if (molecule2.getGlobalBounds().contains(mousey.x,mousey.y)) {MOLECULE = 1;}
                        else if (molecule3.getGlobalBounds().contains(mousey.x,mousey.y)) {MOLECULE = 2;}
                        else if (molecule4.getGlobalBounds().contains(mousey.x,mousey.y)) {MOLECULE = 3;}
                        else if (tracer1.getGlobalBounds().contains(mousey.x,mousey.y)) {TRACER = 0;}
                        else if (tracer2.getGlobalBounds().contains(mousey.x,mousey.y)) {TRACER = 1;}
                        else if (initializeAll.getGlobalBounds().contains(mousey.x, mousey.y) && primed) {
                            TEMPERATURE = sliderT.getSliderValue();
                            NUM = ((int)(sliderN.getSliderValue()/2))*2;
                            DENSITY = sliderD.getSliderValue();
                            STATE = 0;
                            gas = new molec[NUM];
                            initializeMolecules(gas);
                            if (TRACER==0) {tracerX.purge();tracerY.purge();tracerIndex = NUM/2;}
                            t = 0;
                            startAvg = false;
                            eAvg.purge();
                            pAvg.purge();
                            tAvg.purge();}
                    break;}
                case 0: {
                    if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::R) {
                        STATE = -1;
                        primed = false;}
                    break;}
                }
            }}

        // Clear screen
        switch (STATE) {
            case 0: {//simulation
                window.clear();
                generalizedIntegrator(gas);
                for (int i = 0; i < NUM; i++) {
                    ke += 1/2*gas[i].getV()->magsq()*gas[i].getM();
                    ke += 1/2*pow(gas[i].getOmega(),2)*gas[i].getMI();
                    gas[i].disp(&window);}
                ke/=2;
                generalizedCalculation(gas);
                ke = 0;
                pe = 0;
                //showvars
                if (TRACER == 0) {
                    showTracer(&window);
                    if (t %20 ==0) {
                        float ix =gas[tracerIndex].getX();
                        float iy = gas[tracerIndex].getY();
                        tracerX.push(ix*scrnl+effSize*.5);
                        tracerY.push(iy*scrnl+effSize*.5);}}
                if (SHOW_VARS == 0) {
                    sf::RectangleShape background(sf::Vector2f(350,140));
                    background.setPosition(leftalign,upperalign);
                    background.setFillColor(sf::Color::Blue);
                    window.draw(background);
                    energy.setString(estream.str());
                    window.draw(energy);
                    momentum.setString(mstream.str());
                    window.draw(momentum);
                    temp.setString(tstream.str());
                    window.draw(temp);
                    pressure.setString(pstream.str());
                    window.draw(pressure);}
                t++;
                window.display();
                break;}
            case -1: {
                window.clear();
                window.draw(title);
                window.draw(options);
                switch (INTEGRATOR) {
                    case -1: {integrator1.setFillColor(dullWhite);integrator2.setFillColor(dullWhite);break;}
                    case 0: {integrator1.setFillColor(sf::Color::White);integrator2.setFillColor(dullWhite);break;}
                    case 1: {integrator1.setFillColor(dullWhite);integrator2.setFillColor(sf::Color::White);break;}}
                switch (POTENTIAL) {
                    case -1: {potential1.setFillColor(dullWhite);potential2.setFillColor(dullWhite);break;}
                    case 0: {potential1.setFillColor(sf::Color::White);potential2.setFillColor(dullWhite);break;}
                    case 1: {potential1.setFillColor(dullWhite);potential2.setFillColor(sf::Color::White);break;}}
                switch (BOUNDARY) {
                    case -1: {boundary1.setFillColor(dullWhite);boundary2.setFillColor(dullWhite);boundary3.setFillColor(dullWhite);break;}
                    case 0: {boundary1.setFillColor(sf::Color::White);boundary2.setFillColor(dullWhite);boundary3.setFillColor(dullWhite);break;}
                    case 1: {boundary1.setFillColor(dullWhite);boundary2.setFillColor(sf::Color::White);boundary3.setFillColor(dullWhite);break;}
                    case 2: {boundary1.setFillColor(dullWhite);boundary2.setFillColor(dullWhite);boundary3.setFillColor(sf::Color::White);break;}}
                switch (SHOW_VARS) {
                    case -1: {display1.setFillColor(dullWhite);display2.setFillColor(dullWhite);break;}
                    case 0: {display1.setFillColor(sf::Color::White);display2.setFillColor(dullWhite);break;}
                    case 1: {display1.setFillColor(dullWhite);display2.setFillColor(sf::Color::White);break;}}
                switch (SIMULATION) {
                    case -1: {simulation1.setFillColor(dullWhite);simulation2.setFillColor(dullWhite);simulation3.setFillColor(dullWhite);break;}
                    case 0: {simulation1.setFillColor(sf::Color::White);simulation2.setFillColor(dullWhite);simulation3.setFillColor(dullWhite);break;}
                    case 1: {simulation1.setFillColor(dullWhite);simulation2.setFillColor(sf::Color::White);simulation3.setFillColor(dullWhite);break;}
                    case 2: {simulation1.setFillColor(dullWhite);simulation2.setFillColor(dullWhite);simulation3.setFillColor(sf::Color::White);break;}}
                switch (MOLECULE) {
                    case -1: {molecule1.setFillColor(dullWhite);molecule2.setFillColor(dullWhite);molecule3.setFillColor(dullWhite);break;}
                    case 0: {molecule1.setFillColor(sf::Color::White);molecule2.setFillColor(dullWhite);molecule3.setFillColor(dullWhite);molecule4.setFillColor(dullWhite);break;}
                    case 1: {molecule1.setFillColor(dullWhite);molecule2.setFillColor(sf::Color::White);molecule3.setFillColor(dullWhite);molecule4.setFillColor(dullWhite);break;}
                    case 2: {molecule1.setFillColor(dullWhite);molecule2.setFillColor(dullWhite);molecule3.setFillColor(sf::Color::White);molecule4.setFillColor(dullWhite);break;}
                    case 3: {molecule1.setFillColor(dullWhite);molecule2.setFillColor(dullWhite);molecule3.setFillColor(dullWhite);molecule4.setFillColor(sf::Color::White);break;}}
                switch (TRACER) {
                    case -1: {tracer1.setFillColor(dullWhite);tracer2.setFillColor(dullWhite);break;}
                    case 0: {tracer1.setFillColor(sf::Color::White);tracer2.setFillColor(dullWhite);break;}
                    case 1: {tracer1.setFillColor(dullWhite);tracer2.setFillColor(sf::Color::White);break;}}
                window.draw(integrator);
                window.draw(integrator1);
                window.draw(integrator2);
                window.draw(boundary);
                window.draw(boundary1);
                window.draw(boundary2);
                window.draw(boundary3);
                window.draw(potential);
                window.draw(potential1);
                window.draw(potential2);
                window.draw(temperature);
                sliderT.draw(window);
                window.draw(number);
                sliderN.draw(window);
                window.draw(density);
                sliderD.draw(window);
                window.draw(display);
                window.draw(display1);
                window.draw(display2);
                window.draw(simulation);
                window.draw(simulation1);
                window.draw(simulation2);
                window.draw(simulation3);
                window.draw(molecule);
                window.draw(molecule1);
                window.draw(molecule2);
                window.draw(molecule3);
                window.draw(molecule4);
                window.draw(tracer);
                window.draw(tracer1);
                window.draw(tracer2);
                if (POTENTIAL != -1 && INTEGRATOR != -1 && TRACER != -1 && SHOW_VARS != -1 && BOUNDARY != -1 && SIMULATION != -1 && MOLECULE != -1) {primed = true;}
                if (primed) {window.draw(initializeAll);}
                window.display();
                break;}}
        
    }
    return EXIT_SUCCESS;
}
