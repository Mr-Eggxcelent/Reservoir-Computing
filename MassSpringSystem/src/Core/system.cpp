#define _USE_MATH_DEFINES
#include "system.h"
#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <vector>
#include "Eigen/Dense"
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include<map>
#include "../vendor/delaunator-cpp/include/delaunator.hpp"

using namespace Eigen;

//notes
//x=position[0]
//y=position[1]
SpringSystem::SpringSystem(InitialDataValues& data, Camera& camera,unsigned int& Width, unsigned int& Height,bool feedback_state)
    :_data(data),_render(camera, Width, Height),_N(data.N),_feedback_state(feedback_state)
{
    // read out all the data from the data structure
    _num_input_nodes = 0.01 * (_data.input_connectivity_percentage) * _N;
    _num_feedback_nodes = 0.01 * (_data.feedback_connectivity_percentage) * _N;


    //Initialize_Nodes();
    Initialize_Nodes();
    Delaunay_Triangulation_and_Spring_Creation();
    Initialize_Springs();
    Output_For_Plot();

}


void SpringSystem::Initialize_Nodes()
{
    PROFILE_FUNCTION();

    // todo: still needed for debugging
    std::ofstream Node1("src/Output/Node1.csv");
    std::ofstream Node2("src/Output/Node2.csv");
    

    std::ofstream fixed("src/Output/fixednode.txt");
    std::ofstream Initialnodes("src/Output/initial.txt");

    double x;
    double y;

    double x1 = _data.min_x_position;
    double x0 = _data.max_x_position;

    //for fixed nodes.
    int j = 0;
    int k = 0;
    //for bucking center where force is applied
    int c=0;

    std::vector<int>store_line;
    float nodes_that_buckle = 0.3;
    float div = 0;
    float temp = 0;

    for (int i = 0; i < _N; i++)
    {

        if (i < (_N - (_N * nodes_that_buckle)))
        {
            x = Utility::Uniform(_data.min_x_position, _data.max_x_position);

            if (x1 < x)
            {
                x1 = x;
                j = i;

            }

            y = Utility::Uniform(_data.min_y_position, _data.max_y_position);

            if (x0 > x)
            {
                x0 = x;
                k = i;
            }
        }
        else
        {
            div = (_data.max_x_position - _data.min_x_position) / (_N * nodes_that_buckle);
            temp+= _data.min_x_position + div;
            x = temp;

            if (x1 < x)
            {
                x1 = x;
                j = i;

            }

            if (x0 > x)
            {
                x0 = x;
                k = i;
            }

            if ( i== ((_N * (1-nodes_that_buckle)) + ((_N * nodes_that_buckle)/2))-1 )
            {
                c = i;
            }

            store_line.push_back(i);
            y = (_data.max_y_position - _data.min_y_position) / 2;
        }

        Initialnodes << x << "," << y << "\n";
        //The first input_connectivity percent of the nodes are marked as input nodes, and the
        _n.emplace_back(i, x, y, 0);
    }


    fixed << j << "\n";
    fixed << k << "\n";

    //Test to see whether the reason why you're getting those 0 springs is because of fixed nodes.
    _n[j].set_fixed();
    _n[k].set_fixed();

    //buckling functions
    _n[c].set_buckling_center();


    for (int m = 0; m < store_line.size(); m++)
    {
        if (_n[store_line[m]].is_fixed() == false) {
            _n[store_line[m]].set_buckling();
        }
    }


}



void SpringSystem::Delaunay_Triangulation_and_Spring_Creation()
{
    PROFILE_FUNCTION();
 
    double win = 0;
    double wfb = 0;
    int randomnum;


    int num_of_input_nodes = 0;
    int num_of_fixed_nodes = 0;
    int num_of_feedback_nodes = 0;

    std::vector<double> coords;

    for (int i = 0; i < _N; i++)
    {
        //Input weights for the number of input_connectivitiy nodes.
        win = Utility::Uniform(_data.min_input_weight, _data.max_input_weight);
        wfb = Utility::Uniform(_data.min_feedback_weight, _data.max_feedback_weight);

        randomnum = static_cast<int>(Utility::Uniform(0, _N - 1));
   
        //Make sure fixed nodes are not input nodes
        if (i < _num_input_nodes)
        {
            //If it is not a fixed node.
            while (_n[randomnum].is_fixed() || _n[randomnum].is_buckling_node())
            {
                //C++ typecasting rounds down (truncates) but this is fine going from 0 to N-1.
                randomnum = static_cast<int>(Utility::Uniform(0, _N - 1));
            }

            _n[randomnum].init_Input_Node(_data.ux, _data.uy, win);
            num_of_input_nodes++;
        }

        randomnum = static_cast<int>(Utility::Uniform(0, _N - 1));
   
        if(num_of_feedback_nodes < _num_feedback_nodes)
        {
            //If it is not a fixed node.
            while (_n[randomnum].is_fixed() || _n[randomnum].is_buckling_node())
            {
                //C++ typecasting rounds down (truncates) but this is fine going from 0 to N-1.
                randomnum = static_cast<int>(Utility::Uniform(0, _N - 1));
            }

            _n[randomnum].init_Feedback_Node(_data.ux, _data.uy, wfb);
            num_of_feedback_nodes++;

        }

        coords.push_back(_n[i].get_position()[0]);
        coords.push_back(_n[i].get_position()[1]);

    }



    std::cout << "The total number of input nodes is: " << num_of_input_nodes << "\n";
    std::cout << "The total number of feedback nodes is: " << num_of_feedback_nodes << "\n";

    //should not be here will move to when the input nodes are being assigned 
    assign_fb_signal(_data.number_of_signals/_data.order_of_equations);

    Get_Triangles(coords);

}


void SpringSystem::Initialize_Springs()
{
   
    //Spring and damping coefficients
    PROFILE_FUNCTION();

    double k1 = 0;
    double d1 = 0;
    double k3 = 0;
    double d3 = 0;
    double l0 = 0;

    double x0;
    double x1;
    double y0;
    double y1;
    double wout;

    double dist = 0;
    double dist2 = 0;
    double perp_dist_new = 0;
    //maximum possible distance for the initial value as this is a minimisation problem
    double perp_dist_old = _data.max_x_position + _data.max_y_position;

    int connect_node;
    int connect_node2;

    int arraysubscript1 = 0;
    int arraysubscript2 = 0;

    std::ofstream k1output("src/Output/k1output.csv");  // Todo: ofs3, etc. is really bad!
    std::ofstream d1output("src/Output/d1output.csv");
    std::ofstream k3output("src/Output/k3output.csv");
    std::ofstream d3output("src/Output/d3output.csv");
    std::ofstream originallengthoutput("src/Output/originallengthoutput.csv");

    //std::vector<int> node_list;
    //std::vector<int> unconnected_nodes;


    for (unsigned int i = 0; i < _EdgeList.size(); i++)
    {
        //These take the arraysubscripts and disregard the first four points

        arraysubscript1 = _EdgeList[i].first;
        arraysubscript2 = _EdgeList[i].second;

        x0 = _n[arraysubscript1].get_position()[0];
        x1 = _n[arraysubscript2].get_position()[0];
        y0 = _n[arraysubscript1].get_position()[1];
        y1 = _n[arraysubscript2].get_position()[1];


        k1 = Utility::Rand_In_Range_Exp(_data.min_k1, _data.max_k1);
        d1 = Utility::Rand_In_Range_Exp(_data.min_d1, _data.max_d1);
        k3 = Utility::Rand_In_Range_Exp(_data.min_k3, _data.max_k3);
        d3 = Utility::Rand_In_Range_Exp(_data.min_d3, _data.max_d3);


        k1output << k1 << "\n";
        k3output << d1 << "\n";
        d1output << k3 << "\n";
        d3output << d3 << "\n";

        l0 = Utility::Eucl_Dist(x0, y0, x1, y1);
        originallengthoutput << l0 << "\n";

        //Initial value for the output weights. I believe this was never used.
        //wout = Uniform(input_weight_smallest_value, input_weight_largest_value);
        wout = 0;
        _s.emplace_back(arraysubscript1, arraysubscript2, l0, k1, d1, k3, d3, wout);
        

    }

}


void SpringSystem::Reset_Simulation()
{

    double k1_new, d1_new;
    double k3_new;
    double d3_new;

    for (unsigned int i = 0; i < _s.size(); i++)
    {
        k1_new = Utility::Rand_In_Range_Exp(_data.min_k1, _data.max_k1);
        d1_new = Utility::Rand_In_Range_Exp(_data.min_d1, _data.max_d1);
        k3_new = 0;
        d3_new = 0;

        _s[i].reset(k1_new, d1_new, k3_new, d3_new);
        _n[_s[i].get_node_a()].reset_positions();
        _n[_s[i].get_node_b()].reset_positions();

    }

    for (auto& j : _n)
    {
        j.reset_positions();
        j.reset_state();
    }


}



void SpringSystem::calculate_forces(unsigned int& j, double& dt)
{

    _nodea = _s[j].get_node_a();
    _nodeb = _s[j].get_node_b();

    _pos_a = _n[_nodea].get_position();
    _pos_b = _n[_nodeb].get_position();

    _s[j].update2d(_pos_a, _pos_b, dt);

    _F= _s[j].get_force();

    _n[_nodea].apply_force(-_F);
    _n[_nodeb].apply_force(_F);




}

void SpringSystem::assign_fb_signal(int feedback_signal)
{
    //just use a disjoint set this is not the best method
    int split = feedback_signal/_num_feedback_nodes;
    int count = 0;
    int fb_index = 0;


    for (int i = 0; i < _n.size(); i++)
    {
        if (_n[i].is_feedback())
        {
            _assign_signal[i] = fb_index;
            count++;

        }

        if (_n[i].is_feedback() && count == split)
        {
            fb_index++;
            count = 0;
        }

    }

 
}


void SpringSystem::buckle_system(std::vector<Node>::iterator::value_type& l ,double& dt, _STATE& state, bool display)
{

    if (display) {

        switch (state)
        {
        case SpringSystem::_STATE::ORIGINAL: if (l.is_buckling_node() == true) { l.buckle(5.0f, dt); l.set_fixed(); }
                                            _render.render_text("ORIGINAL", 25.0f, 200.0f, 0.25f, glm::vec3(1.0, 1.0f, 1.0f));
                                             break;
        case SpringSystem::_STATE::UP: if (l.is_buckling_node() == true) { l.buckle(7.0f, dt); l.set_fixed(); }
                                       _render.render_text("UP", 25.0f, 200.0f, 0.25f, glm::vec3(1.0, 1.0f, 1.0f));
                                       break;
        case SpringSystem::_STATE::DOWN: if (l.is_buckling_node() == true) { l.buckle(3.0f, dt); l.set_fixed(); }
                                         _render.render_text("DOWN", 25.0f, 200.0f, 0.25f, glm::vec3(1.0, 1.0f, 1.0f));
                                         break;
        }

    }
    else {

        switch (state)
        {
        case SpringSystem::_STATE::ORIGINAL: if (l.is_buckling_node() == true) { l.buckle(5.0f, dt); l.set_fixed(); }
                                             break;
        case SpringSystem::_STATE::UP: if (l.is_buckling_node() == true) { l.buckle(7.0f, dt); l.set_fixed(); }
                                       break;
        case SpringSystem::_STATE::DOWN: if (l.is_buckling_node() == true) { l.buckle(3.0f, dt); l.set_fixed(); }
                                         break;
        }

    }

}

void SpringSystem::update_reset_system(const unsigned int& i,const int& i_2,MatrixXd& input_signal, const MatrixXd& feedback_signal, double& dt, _STATE& state, bool display)
{
    for (auto& l : _n)
    {
        int node_index = 0;

        if (_feedback_state) {
            //Input force to input nodes
            if (l.is_feedback() == true) {
                l.apply_force(Eigen::Vector3d(l.get_w_input() * feedback_signal(i, i_2+_assign_signal[node_index]), 0, 0));
            }
          
        }else {
            //Input force to input nodes assume just single 
            if (l.is_input() == true) {
                l.apply_force(Eigen::Vector3d(l.get_w_input() * input_signal(i,0), 0, 0));
            }
        }

        //Change the node position, velocity and acceleration in response.
        l.update(dt);

        buckle_system(l,dt,state,display);
        
        //At the end of the loop, each node has no force acting on it.
        l.reset_forces();

        node_index++;
    }


}


void SpringSystem::init_draw()
{
    _render.init_text_renderer();

    for (auto& i : _n)
    {
        i.init_render();
    }

    for (auto& i : _s)
    {
        i.init_render();
    }

    _test_line.initBuffer();

}

void SpringSystem::render(glm::mat4& projection, glm::mat4& view)
{
    print_info();

    for (auto& j : _s) {
        j.draw(projection, view, _n);
    }

    for (auto& n : _n)
    {
        std::string input_node = (n.is_input() == true) ? "Input" : " ";
        std::string fb_node = (n.is_feedback() == true) ? "Feedback" : " ";
        std::string node_info = std::to_string(n.get_id()) + " " + input_node+" "+  fb_node;
        _render.render_billboard(node_info, n.get_position()[0], n.get_position()[1] + 0.4f, 0.005, glm::vec3(0.7, 0.8f, 0.3f));
        n.draw(projection, view);
    }


}

void SpringSystem::clearResources()
{
    for (auto& n : _n)
    {
        n.clear_resources();
    }

    for (auto& s : _s)
    {
        s.clear_resources();
    }

}

void SpringSystem::print_info()
{
    std::string details;
    float y_location = 0.25f;

    for (int i = 0; i < _n.size(); i++)
    {
        details = "Node " + std::to_string(_n[i].get_id()) + ": " + std::to_string(_n[i].get_position()[1]);
        _render.render_text(details, 25.0f, y_location, 0.25f, glm::vec3(0.5, 0.8f, 0.2f));
        y_location = y_location + 20;
    }

}


void SpringSystem::Get_Triangles(std::vector<double>&coords)
{
    delaunator::Delaunator d(coords);
    std::map<std::pair<int, int>, int>xy_coords;
    int k = 0;

    for (std::size_t i = 0; i < coords.size(); i+=2)
    {
        xy_coords.insert({ {i, i + 1},k });
        if (i % 2 == 0)
        {
            k++;
        }
    }

    for (std::size_t i = 0; i < d.triangles.size(); i += 3) {
        printf(
            "Triangle points: [[%f, %f], [%f, %f], [%f, %f]]\n",
            d.coords[2 * d.triangles[i]],        //tx0
            d.coords[2 * d.triangles[i] + 1],    //ty0
            d.coords[2 * d.triangles[i + 1]],    //tx1
            d.coords[2 * d.triangles[i + 1] + 1],//ty1
            d.coords[2 * d.triangles[i + 2]],    //tx2
            d.coords[2 * d.triangles[i + 2] + 1] //ty2
        );
    }

    auto find_map = [&](int x_coord,int y_coord) {
        return  xy_coords.find(std::make_pair(x_coord,y_coord))->second;
    };

   
    for (std::size_t i = 0; i < d.triangles.size(); i += 3) {

        _NodeList.first= find_map(2 * d.triangles[i], 2 * d.triangles[i] + 1 );
        _NodeList.second=find_map(2 * d.triangles[i + 1], 2 * d.triangles[i + 1] + 1);
        _EdgeList.push_back(_NodeList);

        _NodeList.first =  find_map(2 * d.triangles[i + 1], 2 * d.triangles[i + 1] + 1);
        _NodeList.second = find_map(2 * d.triangles[i + 2], 2 * d.triangles[i + 2] + 1);
        _EdgeList.push_back(_NodeList);
      

        _NodeList.first =  find_map(2 * d.triangles[i], 2 * d.triangles[i] + 1);
        _NodeList.second = find_map(2 * d.triangles[i + 2], 2 * d.triangles[i + 2] + 1);
        _EdgeList.push_back(_NodeList);
    }

    Utility::Remove_Duplicates(_EdgeList);

}


void SpringSystem::Output_For_Plot()
{
    PROFILE_FUNCTION();

    std::string str = "src/Output/Edges.csv";
    std::ofstream Edges(str);
    std::vector<int>::iterator NodeNums;

    //In Matlab, it does not accept indices of 0 for node graphs.
    //Replace EdgeList with s for consistency.
    for (unsigned int j = 0; j < _s.size() - 1; j++)
    {
        Edges << _s[j].get_node_a() + 1 << "," << _s[j].get_node_b() + 1 << std::endl;
    }

    Edges << _s[_s.size() - 1].get_node_a() + 1 << "," << _s[_s.size() - 1].get_node_b() + 1 << std::endl;

    str = "src/Output/NodeInfo.csv";
    std::ofstream nodesPos(str);

    //just put x,y node position and type together in systems.cpp (still to do)
    unsigned int j = 0;
    while (j < _n.size())
    {
        if (_n[j].is_input()){
            
            nodesPos << _n[j].get_id()+1 <<","<<_n[j].get_position()[0] << "," << _n[j].get_position()[1] << ","<< "I"<< std::endl;

        }else if (_n[j].is_feedback()){
            
            nodesPos << _n[j].get_id()+1 <<","<<_n[j].get_position()[0] << "," << _n[j].get_position()[1] << "," << "FB" << std::endl;

        }else if (_n[j].is_buckling_node()){
            
            nodesPos << _n[j].get_id()+1 <<","<< _n[j].get_position()[0] << "," << _n[j].get_position()[1] << "," << "B" << std::endl;

        }else if (_n[j].is_fixed() && !_n[j].is_buckling_node()){
            
            nodesPos << _n[j].get_id()+1 <<","<< _n[j].get_position()[0] << "," << _n[j].get_position()[1] << "," << "FX" << std::endl;

        }else{
            nodesPos << _n[j].get_id()+1 <<","<< _n[j].get_position()[0] << "," << _n[j].get_position()[1] << "," << "N" << std::endl;
        }

        j++;
    }


}


double SpringSystem::Spring_And_Damping_Coefficient_1(double initial, double finalvalue)
{
    return Utility::Log_10_Uniform(initial, finalvalue);
}


double SpringSystem::Spring_And_Damping_Coefficient_2(double initial, double finalvalue)
{
    return Utility::Uniform(initial, finalvalue);
}

std::vector<Spring> SpringSystem::get_spring_vec()
{
    return _s;
}

std::vector<Node> SpringSystem::get_node_vec()
{
    return _n;
}

unsigned int SpringSystem::Spring_List()
{
    return _EdgeList.size();
}




