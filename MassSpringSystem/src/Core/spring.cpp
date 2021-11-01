//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Original Code by Alan Quille- Bristol University
//https://github.com/AlanQuille
//The code can be found at https://github.com/AlanQuille/Nonlinear-Mass-Spring-System-BR
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "spring.h"
#include <iostream>
#include <cmath>


int Spring::nextID = 0;

Spring::Spring(int node_a, int node_b, double resting_length, double wout)
{
    _id = ++nextID;
    _node_a = node_a;
    _node_b = node_b;
    _resting_length = resting_length;
    _length = resting_length;
    _extension = _length - resting_length;
    _wout = wout;
    _old_extension = _extension;
    
}

Spring::Spring(int node_a, int node_b, double resting_length,double k1, double d1, double wout)
{
    _id = ++nextID;
    _node_a = node_a;
    _node_b = node_b;
    
    _resting_length = resting_length;
    _length = resting_length;
    _extension = _length - resting_length;
    _old_extension = _extension;
  
    _k1 = k1;
    _d1 = d1;
    _wout = wout;
    
}


Spring::Spring(int node_a, int node_b, double resting_length, double k1, double d1, double k3, double d3, double wout)
{
    _id = ++nextID;
    _node_a = node_a;
    _node_b = node_b;
    
    _resting_length = resting_length;
    _length = resting_length;
    _extension = _length - resting_length;
    _old_extension = _extension;
    
    _k1 = k1;
    _d1 = d1;
    _k3 = k3;
    _d3 = d3;
    
    _wout = wout;
  

}

Spring::Spring(const Spring& orig)
{
    _id = orig._id;
    _id = ++nextID;
    _node_a = orig._node_a;
    _node_b = orig._node_b;

    _resting_length = orig._resting_length;
    _length = orig._length;
    _extension = orig._length - orig._resting_length;
    _old_extension = orig._old_extension;

    _k1 = orig._k1;
    _d1 = orig._d1;
    _k3 = orig._k3;
    _d3 = orig._d3;

    _wout = orig._wout;

}

Spring& Spring::operator=(const Spring& orig)
{
    _id = orig._id;
    return(*this);
}



void Spring::update2d(Eigen::Vector3d& node_a_position,Eigen::Vector3d& node_b_position, double dt)
{
    //Calculate nodal displacements in terms of x/y components

    _spring_1 = node_a_position;
    _spring_2 = node_b_position;

    Eigen::Vector3d _displacement = _spring_2 - _spring_1;

    /*
    printf("X1: %f\n", x1);
    printf("X2: %f\n", x2);
    printf("Y1: %f\n", x1);
    printf("Y2: %f\n", x2);
    */

    //Update spring length
    _length = _displacement.norm();

    //Update spring extension and derivative
    _old_extension = _extension;
    _extension = _length - _resting_length;
    _extension_derivative = (_extension - _old_extension) / dt;


    //Use calculate force here to allow option of custom spring force extensions later
    _force_magnitude = calculate_force();
    _force = _force_magnitude * (_displacement / _length);

    /*
    printf("Length: %f\n", length);
    printf("Fx: %f\n", force_x);
    printf("Fy: %f\n", force_y);
    */

    print_debug();

}

double Spring::calculate_force()
{
    double cubic_stiffness = -_k3 * (_extension * _extension * _extension);
    double linear_stiffness = -_k1 * _extension;
    double cubic_damping = -_d3 * (_extension_derivative * _extension_derivative * _extension_derivative);
    double linear_damping = -_d1 * _extension_derivative;


    return linear_stiffness + cubic_stiffness + linear_damping + cubic_damping;

}

//reset 
void Spring::reset(double k1, double d1, double k3, double d3)
{
    //reset instance values
    _length = _resting_length;
    _extension = 0;
    _extension_derivative = 0;
    _force_magnitude = 0;
    _k1 = k1;
    _d1 = d1;
    _k3 = k3;
    _d3 = d3;

}

void Spring::init_render()
{
    _line.initBuffer();
}


void Spring::draw(glm::mat4& projection, glm::mat4& view, std::vector<Node>& nodes)
{
    _line.draw(projection, view, glm::vec3(_spring_1[0], _spring_1[1], 0),
                                 glm::vec3(_spring_2[0], _spring_2[1], 0));

    _line.setColor(glm::vec3(1, 1, 1));
}

void  Spring::clear_resources()
{
    _line.clearResources();

}

void Spring::set_extension(double x1)
{
    _extension = x1;
}

void Spring::set_extension_derivative(double x2)
{
    _extension_derivative = x2;
}

//Return x1 and x2
double Spring::get_length() const
{
    return _length;
}

double Spring::get_resting_length()const
{
    return _resting_length;
}

double Spring::get_extension() const
{
    return _extension;
}

double Spring::get_extension_derivative() const
{
    return _extension_derivative;
}


double Spring::get_force_magnitude() const
{
    return _force_magnitude;
}


//Output the output weight
double Spring::get_Output_Weight() const
{
    return _wout;
}

int Spring::get_id() const
{
    return _id;
}

int Spring::get_node_a() const
{
    return _node_a;
}
int Spring::get_node_b() const
{
    return _node_b;

}

Eigen::Vector3d Spring::get_force() const
{
    return _force;

}


void Spring::set_debug_level(int new_level)
{
    _debug_level = new_level;
}

void Spring::print_debug()
{
    switch (_debug_level)
    {
    case 0:
        break;
    case 1:
        printf("Spring connecting node %d to node %d\n", _node_a, _node_b);
        std::cout << "Length:" << _length<<"\n" ;
        break;
    case 2:
        printf("Spring connecting node %d to node %d\n", _node_a, _node_b);
        printf("Length: %.3f, Force: %.3f\n", _length, _force_magnitude);
        printf("Extension: %.3f, Deriv: %.3f\n", _extension, _extension_derivative);
    }

}

void Spring::print_position()
{
    std::cout << _id << "\t" << _spring_1 << "\t" << _spring_2;
}

void Spring::print_force()
{
    std::cout << _id << "\t" << _force_magnitude << "\t" << _force;
}
void Spring::print_position_and_force()
{
    std::cout << _id << "\t" << _spring_1 <<"\t"<< _spring_2 <<"\t"<< _force_magnitude << "\t" << _force;
}