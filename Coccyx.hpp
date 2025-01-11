//
//  Coccyx.hpp
//  SKELETON
//
//  Created by Owner on 1/6/25.
//

#ifndef Coccyx_hpp
#define Coccyx_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <random>
#include <complex>
#include <iomanip>
using namespace std; 

const long double G= 6.674*pow(10,-11); //N*m^2/kg^2
const long double AU=1.496*pow(10,11); // m
const long double DAYS_IN_SECONDS=86400;
const long double SUNMASS=1.989*pow(10,30); //kg
const long double DT = 86400; //sec (b-a/n) (10800-0/3600) 3 days in seconds



class Body{
private:
    long double mass;
    vector<long double> position;
    vector<long double> velocity;
public:
    Body(long double m, vector<long double> p, vector<long double> v){
        mass=m;
        position=p;
        velocity=v;
    };
    
    long double GetMass(){
        return mass;
    };
    vector<long double> GetPosition(){
        return position;
    }
    vector<long double> GetVelocity(){
        return velocity;
    }
    long double Radius(vector<long double> v,vector<long double>w);
    vector<long double> DistanceVector(vector<long double>v, vector<long double>w);
    vector<long double> GravatationalForce(vector<long double>x, vector<long double>y, long double M, long double m);

    vector<long double> SumOfForces(vector<long double> Rpos,long double m);

    void UpdateVelocityAndPosition(Body b1);
    void PrintResults(vector<long double> r, vector<long double>);
};

bool operator==(Body b1, Body b2);
\
vector<long double> operator+(vector<long double> v, vector<long double> w); //Vector Addition
vector<long double> operator*(long double a, vector<long double> v); //Scalar Multiplication

vector<long double> AUtoMeters(vector<long double> pos); //Convert AU To Meters
vector<long double> MetersToAU(vector<long double> pos); //Convert Meters to AU;
vector<long double> KmPerSecToMetersPerSec(vector<long double> vel); // Km to Meters per sec
vector<long double> MetersPerSecToKmPerSec(vector<long double> vel); // m to Km per sec
long double MAX(long double x, long double y);

#endif /* Coccyx_hpp */



//Excess Code (Please Ignore This Code):
/*
 
 vector<float> SumOfForces(Body b);

 
 //Polymorphism
 class Animal{
 public:
     virtual void MakeASound()=0;
 };

 class Dog:public Animal{
 public:
     void MakeASound(){
         cout<<"The Dog Is Barking";
     }
 };
 class Cat:public Animal{
 public:
     void MakeASound(){
         cout<<"The Cat Is Meowing"<<endl;
     }
 };

 
 struct Branch{
     int data;
     Branch* next=nullptr;
 };

 Branch* AddBrach(int d){
     Branch* result=nullptr;
     result->data=d;
     
     return result;
 }
 
 
 Linked List and Trees:
 
 struct List{
     int data;
     List *next;
 };

 struct Tree{
     int value;
     Tree* Left=nullptr;
     Tree* Right=nullptr;
 };
 List *CreateNode(int d){
     
     List *n=new List();
     n->data=d;
     return n;
 }

 void PrintList(List *n1){
     do{
         cout<<n1->data<<"->";
         n1=n1->next;
     }while(n1 != nullptr);
 }

 Tree *CreateBranch(int d){
     Tree *NewBranch=new Tree();
     NewBranch->value=d;
     
     return NewBranch;
 };
 int main(){
     Tree*r=new Tree();
     Tree*n1=new Tree();
     Tree*n2=new Tree();
     Tree*n3=new Tree();
     Tree*n4=new Tree();


     r->value=5; n1->value=6;
     n2->value=7; n3->value=8; n4->value=9;
     
     r->Left=n1;
     r->Left->Left=n2;
     r->Right=n3;
     r->Right->Left=n4;
 }

 
 
 int size=9;
 PrintShape(size);
 cout<<endl;
 return 0;
 
 //Linked List Code:
 List *root=new List(),*n1=new List(),*n2=new List(),*n3=new List();
 
 root->data=5;
 n1->data=10;
 n2->data=15;
 n3->data=20;
 
 root->next=n1;
 n1->next=n2;
 n2->next=n3;
 n3->next=nullptr;
 
 PrintList(root);
 
 */
