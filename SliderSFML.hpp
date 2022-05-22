//
//  SliderSFML.hpp
//  MD1
//
//  https://github.com/N4G1/Slider-SFML//
//  NOT MINE - TAKEN OFF OF A PUBLIC GIT REPOSITORY AND THEN MODIFIED

#ifndef SliderSFML_hpp
#define SliderSFML_hpp

#include <stdio.h>
#include <SFML/Graphics.hpp>
#include <iostream>
#include <iomanip>
#include <sstream>
class SliderSFML
{
private:
    sf::RectangleShape slider;
    sf::RectangleShape axis;
    sf::Font font;
    sf::Text text;
    std::stringstream val, hi, lo;
    float minValue;
    float maxValue;
    int xCord;
    int yCord;
    int axisWidth;
    int axisHeight;
    int sliderWidth;
    int sliderHeight;
    float sliderValue;
    bool datType; //true is int, false is float
public:
    SliderSFML(int x, int y);
    sf::Text returnText(int x, int y, std::string z, int fontSize);
    void create(float min, float max, bool isInt);
    void logic(sf::RenderWindow &window);
    float getSliderValue();
    void setSliderValue(float newValue);
    void setSliderPercentValue(float newPercentValue);
    void draw(sf::RenderWindow & window);
};
#endif /* SliderSFML_hpp */
