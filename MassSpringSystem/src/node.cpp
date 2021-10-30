//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Original Code by Alan Quille- Bristol University
//https://github.com/AlanQuille
//The code can be found at https://github.com/AlanQuille/Nonlinear-Mass-Spring-System-BR
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "node.h"


Node::Node(int id,double x_position, double y_position, double z_position, double mass)
{
    _id = id;

    _position << x_position, y_position, z_position;

    _velocity = Eigen::Vector3d::Zero();

    _initial_position << x_position, y_position, z_position;

    _mass = mass;

}

Node::Node(const Node& orig) {

    _id = orig._id;

    _position = orig._position;

    _velocity = orig._velocity;

    _initial_position = orig._initial_position;

    _mass = orig._mass;
}

Node& Node::operator=(const Node& orig) {
    _id = orig._id;
    return(*this);
}

//If this function is activated than the node is input node
void Node::init_Input_Node(double ux, double uy, double win)
{
    _input_force << ux, uy, 0;
    _input_node = true;
    _w_input = win;
    std::cout << "win is: " << _w_input<<"\n";
}

void Node::init_Feedback_Node(double ux, double uy, double wfb)
{
    _feedback_force << ux, uy, 0;
    _feedback_node = true;
    _w_feedback = wfb;
    std::cout << "wfb is: " << _w_feedback<<"\n";
}


//This is the function that incrementally changes the Node position in the next timestep;
void Node::apply_force(const Eigen::Vector3d& F)
{
    if (_fixed_node == false)
    {
        _input_force += F;
    }
}

//Thanks to Jorge Rodriguez tutorial  https://www.youtube.com/watch?v=qJq7I2DLGzI
double Node::approach(double flGoal, double flCurrent, double dt)
{
    //temporary empty
    double difference = flGoal - flCurrent;

    if (difference > dt)
    {
        return flCurrent += dt;
    }

    if (difference < dt)
    {
        return flCurrent -= dt;
    }

    return flGoal;

}

void Node::buckle(double target_pos, double dt)
{
    double delta_y = target_pos - get_position()[1];
    double goal_dist = sqrt(delta_y * delta_y);

    if (goal_dist >0.25) {

        _velocity[1] = approach(delta_y, _velocity[1], dt * 100);
        _position[1] = _position[1] + dt * _velocity[1];
        _velocity[1] = _velocity[1] + dt;

    }
   
}


void Node::update(double dt)
{
    if (_fixed_node == false)
    {
        //Standard Euler integration - maybe add Runge Kutta later
        _acceleration = _input_force / _mass;
        
        _velocity += dt * _acceleration;

        _position += dt * _velocity;

    }

}


void Node::reset_positions()
{
    _position = _initial_position;

}

void Node::reset_forces()
{
    _input_force = Eigen::Vector3d::Zero();
}


void Node::reset()
{
    reset_state();
    reset_forces();
    reset_positions();
}

void Node::reset_state()
{
    _velocity = Eigen::Vector3d::Zero();
    _acceleration = Eigen::Vector3d::Zero();

}

void Node::set_fixed()
{
    _fixed_node = true;
}

void Node::set_buckling()
{
    _buckling_node = true;
}

void Node::set_buckling_center()
{
    _central_buckling_node = true;
}

void Node::set_mass(double new_mass)
{
    _mass = new_mass;
}

void Node::init_render()
{
    _sphere.initBuffer();
}


void Node::draw(glm::mat4& projection, glm::mat4& view)
{
    if (is_fixed() == true)
    {
        _sphere.draw(projection, view, glm::vec3(get_position()[0], get_position()[1], get_position()[2]), glm::vec4(1.0f, 0.0f, 0.0f, 1.0f));
    }
    else if (this->is_input() == true )
    {
        _sphere.draw(projection, view, glm::vec3(get_position()[0], get_position()[1], get_position()[2]), glm::vec4(0.0f, 1.0f, 0.0f, 1.0f));
    }
    else
    {
        _sphere.draw(projection, view, glm::vec3(get_position()[0], get_position()[1], get_position()[2]), glm::vec4(0.0f, 0.0f, 1.0f, 1.0f));
    }


}


void Node::clear_resources()
{
    _sphere.clearResources();
}


bool Node::is_fixed()
{
    return _fixed_node;
}

bool Node::is_buckling_node()
{
    return _buckling_node;
}


bool Node::is_central_buckling()
{
    return _central_buckling_node;
}


double Node::get_mass() const
{
    return _mass;
}

int Node::get_id() const
{
    return _id;
}

double Node::get_w_input() const
{
    return _w_input;
}

double Node::get_w_feedback() const
{
    return _w_feedback;
}


bool Node::is_input()
{
    return _input_node;
}

bool Node::is_feedback()
{
    return _feedback_node;
}

Eigen::Vector3d Node::get_position()
{
    return _position;
}

void Node::print_position()
{
    std::cout << _id << "\t" << _position;
}

void Node::print_positions_and_velocities()
{
    std::cout << _id << "\t" << _position<<"\t"<<_velocity;
}


void Node::change_update_check()
{
    _update_check = 1;
}

bool Node::return_update_check()
{
    return _update_check;
}

void Node::zero_update_check()
{
    _update_check = 0;
}


