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

    //Earths Motion Implementation:
    
    Body b(5.972*pow(10,24),{1.0,0,0},{0.0,29.78,0.0});
    b.UpdateVelocityAndPosition(b);
    return 0;
}







