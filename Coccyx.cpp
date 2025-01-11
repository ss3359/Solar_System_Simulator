//
//  Coccyx.cpp
//  SKELETON
//
//  Created by Owner on 1/6/25.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <random>
#include <complex>
#include <algorithm>
#include <iomanip>
#include "Coccyx.hpp"


using namespace std;
Body Mercury(3.301e23,{0.307,0,0},{0.0,47.87,0.0});
Body Venus(4.867e24,{0.718,0,0},{0.0,35.02,0.0});
Body Earth(5.972e24,{1.0,0,0},{0.0,29.78,0.0});
Body Mars(6.417e24,{1.381,0,0},{0.0,24.07,0.0});
Body Jupiter(1.898e27,{4.95,0,0},{0.0,13.07,0.0});
Body Saturn(5.683e26,{9.0,0,0},{0.0,9.69,0.0});
Body Uranus(8.681e25,{19.0,0,0},{0.0,6.80,0.0});
Body Neptune(1.024e26,{30.0,0,0},{0.0,5.43,0.0});
vector<Body> Planets={Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune};

long double Body:: Radius(vector<long double> v,vector<long double>w){
    long double result=0;
    for(int i=0; i<v.size(); i++){
        result+=pow(w[i]-v[i],2);
    }
    return sqrt(MAX(result,1e-10));
}


vector<long double> Body:: DistanceVector(vector<long double>v, vector<long double>w){
    vector<long double> result;
    for(int i=0; i<v.size(); i++){
        result.push_back(w[i]-v[i]);
    }
    return result;
}

vector<long double> Body::GravatationalForce(vector<long double>x, vector<long double>y, long double M, long double m){
    vector<long double> distance_vector(3,0);
    long double r=0.0f;
    
    x=AUtoMeters(x); y=AUtoMeters(y);
    
    for(int i=0; i<x.size(); i++){
        distance_vector[i]=y[i]-x[i];
        r+=pow(distance_vector[i],2);
    }
    if(r<1e-6){
        return vector<long double>(3,0);
    }
    vector<long double>unit_vector(3,0);
    for(int i=0; i<distance_vector.size(); i++){
        unit_vector[i]=distance_vector[i]/r;
    }

    long double ForceScalar= G*(M*m)/(r*r);
    for(int i=0; i<unit_vector.size(); i++){
        unit_vector[i]*=ForceScalar;
}
    cout<<"Gravitational Force:";
    for(long double g: unit_vector)cout<<g<<" ";
    cout<<endl;
        
    
    long double force_mag=0.0f;
    for(int i=0; i<unit_vector.size(); i++){
        force_mag+=unit_vector[i]*unit_vector[i];
    }
    force_mag=sqrt(force_mag);
    cout<<"Force Magnitude: "<<force_mag<<endl;
        return unit_vector;
}




/*
 The differential equation for the RK Method is the Following
 dv/dt = F(t)/m
 */

/*
 This function represents the sum of the gravatational forces acting on the single celestial body, b.
 */
bool operator==(Body b1, Body b2){
    long double tol = 1e-6;
    return (abs(b1.GetMass()-b2.GetMass())<tol)
    &&(b1.GetPosition()==b2.GetPosition())
    &&(b1.GetVelocity()==b2.GetVelocity());
}

vector<long double> operator+(vector<long double> v, vector<long double> w){
    vector<long double> result;
    
    for(int i=0; i<v.size(); i++){
        result.push_back(v[i]+w[i]);
    }
    return result;
}

vector<long double> operator*(long double a, vector<long double> v){
    vector<long double> result;
    
    for(int i=0; i<v.size(); i++){
        result.push_back(a*v[i]);
    }
    return result;
}

vector<long double> AUtoMeters(vector<long double> pos){
    return AU*pos;
} //Convert AU To Meters
vector<long double> MetersToAU(vector<long double> pos){
    return(1.0/AU)*pos;
} //Convert Meters to AU;
vector<long double> KmPerSecToMetersPerSec(vector<long double> vel){
    return (1000.0)*vel;
}// Km to Meters per sec
vector<long double> MetersPerSecToKmPerSec(vector<long double> vel){
    return (1.0/1000.0)*vel;
} // m to Km per sec
long double MAX(long double x, long double y){
    if(x>=y)
        return x;
    else
        return y;
}




vector<long double> Body::SumOfForces(vector<long double> Rpos,long double m){
    if(m==0){
        cout<<"Error Mass is Zero!";
        return vector<long double>(3,0);
    }
    vector<long double> F(3,0);
    for(Body p: Planets){
        if(p==*this)
                continue;
        vector<long double> G=GravatationalForce(Rpos, p.GetPosition(), m, p.GetMass());

        for(int i=0; i<F.size();i++){
            F[i]+=G[i];
        }
    }
    cout<<"Force Vector: ";
    for(long double f: F) cout<<f<<" ";
    cout<<endl;
    return F;
}

void Body::UpdateVelocityAndPosition(Body b1){
    vector<long double> Rpos=AUtoMeters(b1.GetPosition()); //Position Vector;
    vector<long double> v=KmPerSecToMetersPerSec(b1.GetVelocity()); //Velocity Vector;
    long double m=b1.GetMass();
    
    
    for(Body b: Planets){
        if(b==*this) continue;
        vector<long double> d = AUtoMeters(DistanceVector(b.GetPosition(), this->GetPosition()));
        
        cout<<endl;
        
        cout<<"Distance Vector:";
        for(long double i: d) cout<<i<<"  ";
        cout<<endl;
        
        //Constants For Runge-Kutta Method
        for(int n=0; n<1; n++){
            
            vector<long double> k1=v;
            vector<long double> l1=SumOfForces(Rpos, m);
            
            vector<long double> k2= v+(DT/2.0)*l1;
            vector<long double> Rpos_k2=Rpos+(DT/2.0)*k1;
            vector<long double> l2=SumOfForces(Rpos_k2, m);

            vector<long double> k3=v+(DT/2.0)*l2;
            vector<long double> Rpos_k3=Rpos+(DT/2.0)*k2;
            vector<long double> l3=SumOfForces(Rpos_k3, m);
            
            vector<long double> k4= v+(DT)*l3;
            vector<long double> Rpos_k4=Rpos+DT*k3;
            vector<long double> l4=SumOfForces(Rpos_k4,m);
    
            cout<<"l1: ";
            for(long double f: l1) cout<<f<<" ";
            cout<<endl;
            
            cout<<"l2: ";
            for(long double f: l2) cout<<f<<" ";
            cout<<endl;
            
            cout<<"l3: ";
            for(long double f: l2) cout<<f<<" ";
            cout<<endl;
            
            cout<<"l4: ";
            for(long double f: l2) cout<<f<<" ";
            cout<<endl;
            
            for(int i=0; i<3; i++){
                v[i]=v[i]+(DT/6.0)*(l1[i]+(2.0*l2[i])+(2.0*l3[i])+l4[i]);
                Rpos[i]=Rpos[i]+(DT/6.0)*(k1[i]+(2.0*k2[i])+(2.0*k3[i])+k4[i]);
            }
            cout<<endl;
            PrintResults(Rpos,v);
        }
        
    }
}
void Body::PrintResults(vector<long double> r, vector<long double> v){
    cout<<setprecision(8);
    cout<<"Position Vector: (";
    for(int i=0; i<r.size(); i++){
        if(i==r.size()-1)
            cout<<r[i];
        else
            cout<<r[i]<<",";
    }
    cout<<")"<<endl;
    
    cout<<"Veclocity Vector: (";
    for(int i=0; i<v.size(); i++){
        if(i==v.size()-1)
            cout<<v[i];
        else
            cout<<v[i]<<",";
    }
    cout<<")"<<endl;
}

