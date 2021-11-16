#pragma once

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "spring.h"
#include "Instrumentor.h"
#include "utility.h"
#include "../Graphics/render_text.h"

using namespace Eigen;


struct InitialDataValues
{
    int     N;   // number of Nodes
    double ux;   // First input values in x direction  TODO: Really needed?
    double uy;   // First input values in y direction
    //double uz;   // First input values in z direction

    int order_of_equations;
    int number_of_signals;
    int number_of_equations;

    double input_connectivity_percentage;  // [0,1] percentage of nodes that receive input
    double feedback_connectivity_percentage;  // [0,1] percentage of nodes that receive input

    // lower and upper range for input weights
    double min_input_weight;
    double max_input_weight;

    // lower and upper range for input weights
    double min_feedback_weight;
    double max_feedback_weight;

    // Parameters to set area where nodes can be placed in
    // Todo: maybe change to min_x... and max... also below with probabilties
    double min_x_position;
    double max_x_position;

    double min_y_position;
    double max_y_position;

    //double min_z_position;
    //double max_z_position;

    double t0;      // time for first time step [s]
    double tmax;    // maximum time step [s]
    double dt;      // time step in seconds

    //The min and max values for k1, k3, d1, d3 values, either log uniform or uniform depending.
    double min_k3;
    double max_k3;
    double min_d3;
    double max_d3;

    double min_k1;
    double max_k1;
    double min_d1;
    double max_d1;
};






class SpringSystem
{


private:
    InitialDataValues _data;
    unsigned int _N;          // Number of mass points
    double _num_input_nodes;    // Number of input nodes
    double _num_feedback_nodes; // Number of feedback nodes

    std::vector<Spring> _s;      // List of all springs
    std::vector<std::pair<double,double>> _EdgeList;   
    std::pair<double,double> _NodeList;     
    std::vector<double> _EdgeNodeList;

    //use a modified disjoint set it will work better than a map
    std::map<int, int>_assign_signal;
    bool _feedback_state;

    RenderText _render;
    LineRender _test_line;

    bool _input_node;
    bool _feedback_node;

    unsigned int _nodea = 0;
    unsigned int _nodeb = 0;

    Eigen::Vector3d _pos_a;
    Eigen::Vector3d _pos_b;

    double _vector_x = 0;
    double _vector_y = 0;

    Eigen::Vector3d _F;
    double _Fx = 0;
    double _Fy = 0;

  

public:
    std::vector<Node> _n;        // List of all nodes
    //enum class _STATE {ORIGINAL=0,DOWN=1,UP=2,DEFAULT=3};
    enum class _STATE { ORIGINAL=0, UP=1, DOWN=2, DEFAULT=3 };
     //enum class _STATE { UP=0, ORIGINAL=1, DOWN=2, DEFAULT=3 };
    _STATE _last_state;

public:

    SpringSystem(InitialDataValues& data,Camera& camera,unsigned int& Width,unsigned int& Height, bool feedback_state);

    //This initializes the nodes and puts in appropriate values for the ranges and the weightsstd
   // Todo: change variable names here
    void Initialize_Nodes();

    //This does the delaunay triangulation for the two dimensional case and creates the springs for the reservoir computer, not the radial spider web
    void Delaunay_Triangulation_and_Spring_Creation();

    //Create springs for reservoir computer nonlinear mass spring system
    void Initialize_Springs();

    //This changes position of springs and nodes dynamically in time.
    void Reset_Simulation();


    //update
    void calculate_forces(unsigned int& j, double& dt);

    void assign_fb_signal(int feedback_signal);

    void buckle_system(std::vector<Node>::iterator::value_type& l, double& dt, _STATE& state, bool display);

    //system
    /*template <typename Derived>
    const EigenBase<Derived>& feedback_signal*/
    void update_reset_system(const unsigned int& i,const int& i_2, MatrixXd& input_signal, const MatrixXd& feedback_signal, double& dt, _STATE& state, bool display = false);

    //Init the buffers before drawing
    void init_draw();

    //draw to screen
    void render(glm::mat4& projection, glm::mat4& view);

    //Init the buffers before drawing
    void clearResources();

    //Debug to screen
    void print_info();

    //After delaunay triangulation, how many connecting edges hence springs.
    // Todo: Should the ouput be really double and not unsigned int?
    unsigned int Output_No_of_Edges();

    //This outputs the spring and node positions for the reservoir computer
    void Output_Spring_And_Node_Positions();

    //Get the triangles from the Delaunay Triangulation.
    void Get_Triangles(std::vector<double>&coords);

    void Output_For_Plot(std::string system_pos);

    //Expression of log 10 uniform for spring and damping coefficients
    double Spring_And_Damping_Coefficient_1(double initial, double finalvalue);
    double Spring_And_Damping_Coefficient_2(double initial, double finalvalue);

    std::vector<Spring>get_spring_vec();

    //Return node vector
    std::vector<Node> get_node_vec();

    //Return number of edges from the triangle.
    unsigned int Spring_List();

 


};
