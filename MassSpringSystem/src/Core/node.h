//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Original Code by Alan Quille- Bristol University
//https://github.com/AlanQuille
//The code can be found at https://github.com/AlanQuille/Nonlinear-Mass-Spring-System-BR
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma once
#include <iostream>
#include "../Graphics/shapes.h"
#include "Eigen/Dense"

//notes
//x=position[0]
//y=position[1]
class Node
{
	int _id;
	double _mass = 1;

	//Original positions
	Eigen::Vector3d _initial_position;
	// Cartesian coordinates
	Eigen::Vector3d _position;
	// Velocity
	Eigen::Vector3d _velocity;
	// Acceleration 
	Eigen::Vector3d _acceleration = Eigen::Vector3d::Zero();


	// Force applied to the Node
	Eigen::Vector3d _force = Eigen::Vector3d::Zero();
	// Spring force applied to the Node
	Eigen::Vector3d _spring_force = Eigen::Vector3d::Zero();
	// Input force applied to the Node
	Eigen::Vector3d _input_force = Eigen::Vector3d::Zero();
	// Feedback force applied to the Node
	Eigen::Vector3d _feedback_force = Eigen::Vector3d::Zero();
	// Buckling force applied to the Node
	Eigen::Vector3d _buckling_force = Eigen::Vector3d::Zero();

	// The velocity of the nodes initially it is zero.
	// int connections;  // TODO: Still needed
	bool _input_node = false;
	bool _feedback_node = false;
	bool _buckling_node = false;
	bool _central_buckling_node = false;


	// Node is globally fixed or not
	bool _fixed_node = false;
	double _w_input = 0;
	double _w_feedback = 0;

	bool _moving = false;
	//A bool check to determine whether a node has been updated or not.
	bool _update_check = 0;
	
	SphereRender _sphere;
	

public:

	//Default constructor, with caretesian coordinates
	Node(int id,double x_position, double y_position, double z_position, double mass = 1.0);
	Node(const Node& orig);
	Node& operator=(const Node& orig);

	// Function to convert Node to node that receives and an input force
	void init_Input_Node(double ux, double uy, double win);
	// Function to convert Node to node that receives and an feedback force
	void init_Feedback_Node(double ux, double uy, double wfb);

	//This is the function that incrementally changes the nodes position in the next timestep;
	// maybe instead of "change" use "update"
	void apply_spring_force(const Eigen::Vector3d& F);
	void apply_input_force(const Eigen::Vector3d& F);
	void apply_feedback_force(const Eigen::Vector3d& F);
	void apply_buckling_force(const Eigen::Vector3d& F);

	void buckle(double target_pos, double dt, bool& key_lock);
	void update(double dt);

	//Save original position in x and y of nodes.
	void reset_positions();
	void reset_forces();
	void reset();
	void reset_state();


	// Make node globally fixed
	void set_fixed();
	void set_buckling();
	void set_buckling_center();
	void set_mass(double new_mass);


	//Draw Nodes
	void init_render();
	void draw(glm::mat4& projection, glm::mat4& view);
	void clear_resources();


	//Returns fixed node
	bool is_fixed();
	bool is_buckling_node();
	bool is_central_buckling();

	double get_mass()const;
	int get_id()const;
	double get_w_input() const;
	double get_w_feedback() const;
	
	// Checks if node is input node
	bool is_input();
	bool is_feedback();

	Eigen::Vector3d get_position();

	//Show current position
	void print_position();
	void print_positions_and_velocities();

	//At every timestep, a node should be changed only once. A node should not be changed at every tiem
	void change_update_check();

	//Check to see if the node has been modified or not
	bool return_update_check();

	//Check to see if the node has been modified or not
	void zero_update_check();


};

