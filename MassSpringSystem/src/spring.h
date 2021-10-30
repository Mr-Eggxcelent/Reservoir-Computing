//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Original Code by Alan Quille- Bristol University
//https://github.com/AlanQuille
//The code can be found at https://github.com/AlanQuille/Nonlinear-Mass-Spring-System-BR
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma once
#include<iostream>
#include "../core/shapes.h"
#include "node.h"

const double K1_DEFAULT = 5;
const double K3_DEFAULT = 0;
const double D1_DEFAULT = .1;
const double D3_DEFAULT = 0;

class Spring
{

private:
    int _id;
    static int nextID;

    double _resting_length;  // resting length l0
    double _length;   // current length.
    double _extension = 0; // difference between current length and resting length x1
    double _extension_derivative = 0; // derivative of extension x2
    double _old_extension = 0;


    // spring constants
    // force from the spring stiffness based on x1 is
    // p = x1^3*k3 + x1*k1
    double _k1=K1_DEFAULT;
    double _k3=K3_DEFAULT;

    // damping constants
    // force from the damping based on x2 (change of length over time)
    // q = x2^3*d3 + x3*d1
    double _d1=D1_DEFAULT;
    double _d3= D3_DEFAULT;

    // Two nodes that are connected by this spring
    int _node_a;
    int _node_b;

    //output weight assigned to this spring
    double _wout;

    //Total force applied by this spring
    double _force_magnitude = 0;
    //Force components
    Eigen::Vector3d _force = Eigen::Vector3d::Zero();
    int _debug_level = 0;

    //These represent the end points of the spring; store to make printing consistent with nodes
    Eigen::Vector3d _spring_1 = Eigen::Vector3d::Zero();
    Eigen::Vector3d _spring_2 = Eigen::Vector3d::Zero();

    LineRender _line;
   

public:
    // Default constructor to load in spring and damping coefficeints
   // Todo: We should make this classe more general - the stiffness and damping functiosn should be overloaded
   // So we can implement variations of that
    Spring(int node_a, int node_b, double resting_length, double wout); //Make a spring with default stiffness params
    Spring(int node_a, int node_b, double resting_length, double k1, double d1, double wout); //Make a linear spring
    Spring(int node_a, int node_b, double resting_length, double k1, double d1, double k3, double d3, double wout);
    
    Spring(const Spring& orig);
    Spring& operator=(const Spring& orig);


    void update2d(Eigen::Vector3d& node_a_position, Eigen::Vector3d& node_b_position, double dt);

    double calculate_force();

    //reset 
    void reset(double k1, double d1, double k3, double d3);

    void init_render();
    void draw(glm::mat4& projection, glm::mat4& view, std::vector<Node>& nodes);
    void clear_resources();

    void set_extension(double x1);
    void set_extension_derivative(double x2);
    
    //Return x1 and x2
    double get_length() const;
    double get_resting_length()const;
    double get_extension() const;
    double get_extension_derivative() const;

    double get_force_magnitude() const;   

    //Output the output weight
    double get_Output_Weight()const;
    int get_id()const;
    int get_node_a()const;
    int get_node_b()const;

    Eigen::Vector3d get_force() const;

    void set_debug_level(int new_level);
    void print_debug();

    void print_position();
    void print_force();
    void print_position_and_force();

};