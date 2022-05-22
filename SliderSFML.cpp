//
//  SliderSFML.cpp
//  MD1
//
//  Created by Connor Blake on 11/25/21.
//  Copyright © 2021 Connor Blake. All rights reserved.
//

#include "SliderSFML.hpp"
SliderSFML::SliderSFML(int x, int y)
{
    xCord = x;
    yCord = y;
    axisHeight = 10;
    axisWidth = 400;
    sliderWidth = 20;
    sliderHeight = 30;
    
    if (!font.loadFromFile("/Users/connorblake/Documents/SORTOld/RampCurves/RampCurves/arial.ttf"))
        std::cout << "Error loading font\n";

    text.setFont(font);
    text.setFillColor(sf::Color::White);

    axis.setPosition(x, y);
    axis.setOrigin(0, axisHeight / 2);
    axis.setSize(sf::Vector2f(axisWidth, axisHeight));
    axis.setFillColor(sf::Color(63,63,63));
    slider.setPosition(x, y);
    slider.setOrigin(sliderWidth / 2, sliderHeight / 2);
    slider.setSize(sf::Vector2f(sliderWidth, sliderHeight));
    slider.setFillColor(sf::Color(192,192,192));
}

sf::Text SliderSFML::returnText(int x, int y, std::string z, int fontSize)
{
    text.setCharacterSize(fontSize);
    text.setPosition(x, y);
    text.setString(z);
    return text;
}

void SliderSFML::create(float min, float max, bool isInt)
{
    sliderValue = min+(max-min)*3/8;
    minValue = min;
    maxValue = max;
    datType = isInt;
    lo << std::fixed << std::setprecision(2) << minValue;
    hi << std::fixed << std::setprecision(2) << maxValue;
}

void SliderSFML::logic(sf::RenderWindow &window)
{
    if (slider.getGlobalBounds().contains(sf::Mouse::getPosition(window).x, sf::Mouse::getPosition(window).y)
        && sf::Mouse::isButtonPressed(sf::Mouse::Button::Left))
    {
        if (sf::Mouse::getPosition(window).x >= xCord && sf::Mouse::getPosition(window).x <= xCord + axisWidth)
        {
            slider.setPosition(sf::Mouse::getPosition(window).x, yCord);
            sliderValue = (minValue + ((slider.getPosition().x - xCord) / axisWidth * (maxValue - minValue)));
        }
    }
}

float SliderSFML::getSliderValue()
{
    return sliderValue;
}

void SliderSFML::setSliderValue(float newValue)
{
    if (newValue >= minValue && newValue <= maxValue)
    {
        sliderValue = newValue;
        float diff = maxValue - minValue;
        float diff2 = newValue - minValue;
        float zzz = axisWidth / diff;
        float posX = zzz*diff2;
        posX += xCord;
        slider.setPosition(posX, yCord);
    }
}

void SliderSFML::setSliderPercentValue(float newPercentValue)
{
    if (newPercentValue >= 0 && newPercentValue <= 100)
    {
        sliderValue = newPercentValue / 100 * maxValue;
        slider.setPosition(xCord + (axisWidth*newPercentValue / 100), yCord);
    }
}

void SliderSFML::draw(sf::RenderWindow &window)
{
    logic(window);
    if (datType) {
        window.draw(returnText(xCord - 10, yCord + 5, std::to_string((int)minValue), 20));
        window.draw(returnText(xCord + axisWidth - 10, yCord + 5, std::to_string((int)maxValue), 20));
        window.draw(returnText(slider.getPosition().x - sliderWidth, slider.getPosition().y - sliderHeight,
            std::to_string((int)sliderValue), 15));
    }
    else{
        window.draw(returnText(xCord - 10, yCord + 5, lo.str(), 20));
        window.draw(returnText(xCord + axisWidth - 10, yCord + 5, hi.str(), 20));
        val.str(std::string());
        val.clear();
        val << std::fixed << std::setprecision(2) << sliderValue;
        window.draw(returnText(slider.getPosition().x - sliderWidth, slider.getPosition().y - sliderHeight,val.str(), 15));}
    window.draw(axis);
    window.draw(slider);

}
