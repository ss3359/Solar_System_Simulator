//
//  main.cpp
//  SKELETON
//
//  Created by Owner on 11/18/24.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <random>
#include <complex>
#include <iomanip>

#include "Coccyx.hpp"

using namespace std;


/*
I am going to code a Solar System Simulation with script first. Then, I will use an 3d graphics library to simulate a Solar System. First, I need to define the initial conditions for the respective planets of the solar system. The intial conditions are for the mass, position, and velocity of the celestial object.
 */

//Global Variables

//Celestial Bodies
/*
 The units for the initial conditions for the planets are as follows
 Mass(kg), Position(x,y,z)(AUs), Velocity (vx, vy,vz)(km/s)
 */


int main(){

    Body Mercury(3.301e23,{0.307,0,0},{0.0,47.87,0.0});
    Body Venus(4.867e24,{0.718,0,0},{0.0,35.02,0.0});
    Body Earth(5.972e24,{1.0,0,0},{0.0,29.78,0.0});
    Body Mars(6.417e24,{1.381,0,0},{0.0,24.07,0.0});
    Body Jupiter(1.898e27,{4.95,0,0},{0.0,13.07,0.0});
    Body Saturn(5.683e26,{9.0,0,0},{0.0,9.69,0.0});
    Body Uranus(8.681e25,{19.0,0,0},{0.0,6.80,0.0});
    Body Neptune(1.024e26,{30.0,0,0},{0.0,5.43,0.0});
    Body Planets[8]={Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune};

    
    
    //Earths Motion Implementation:
    
    Body b(Earth.GetMass(),Earth.GetPosition(),Earth.GetVelocity());
    b.UpdateVelocityAndPosition(b,Planets);
    return 0;
}







