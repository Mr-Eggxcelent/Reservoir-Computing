#define _USE_MATH_DEFINES

#include "simulation.h"
#include <intrin.h>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <Eigen/QR>
#include "filereader.h"
#include<memory>
#include "../Graphics/camera.h"

//graphing
#include "../Graphics/gnuplot-iostream.h"


using namespace Eigen;

void no_feedback_generator();
void feedback_generator();
void readParameter(std::map <std::string, double>&,std::string);


#ifdef _WIN32
unsigned long long rdtsc() {
    return __rdtsc();
}

//  Linux/GCC
#else
unsigned long long rdtsc() {
    unsigned int lo, hi;
    __asm__ __volatile__("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}
#endif





int main(int argc, char** argv)
{

    auto begin = std::chrono::high_resolution_clock::now();

    no_feedback_generator();
   //feedback_generator();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "The time it took for the programme to run in total in milliseconds: ";
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms \n";


}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////Volterra and NARMA////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void no_feedback_generator()
{

    InitialDataValues data;
    srand(rdtsc());

    std::vector<double>volterra;
    std::vector<double>second_order;
    std::vector<double>narma;
    std::vector<double>input;
    
    std::cout << "-- Start ---------------------------------------- \n";

    auto read_lambda = [](std::string file_name, std::vector<double>& vec_signal,int column = 1) {
        FileReader input_signal(file_name, vec_signal);
        input_signal.file_read(column);
    };

    std::thread input_thread(read_lambda, "src/Data/inputsignal.csv", std::ref(input),1);
    std::thread volterra_thread(read_lambda, "src/Data/volterra.csv", std::ref(volterra),2);
    std::thread second_order_thread(read_lambda, "src/Data/2ndOrder.csv", std::ref(second_order),0);
    std::thread NARMA_thread(read_lambda, "src/Data/NARMA.csv", std::ref(narma),1);

    volterra_thread.join();
    second_order_thread.join();
    NARMA_thread.join();
    input_thread.join();
    
    MatrixXd Target_Signals(second_order.size(), 3);
    MatrixXd Input_Signals(input.size(), 1);
    
  
    Target_Signals.col(0)  = Eigen::Map<Eigen::VectorXd>(volterra.data(), volterra.size());
    Target_Signals.col(1)  = Eigen::Map<Eigen::VectorXd>(second_order.data(), second_order.size());
    Target_Signals.col(2)  = Eigen::Map<Eigen::VectorXd>(narma.data(), narma.size());
    Input_Signals .col(0)  = Eigen::Map<Eigen::VectorXd>(input.data(), input.size());

    std::map <std::string, double>var_map;
    readParameter(var_map, "src/Data/DataNoFeedback.txt");

    double wash_out_time = var_map["wash_out_time"];
    double learning_time = var_map["learning_time"];
    double learning_time_test = var_map["learning_time_test"];

    // Setting parameters for simulation
    // This should be possible to read in from a text file
    data.N = static_cast<int>(var_map["N"]);
    data.ux = var_map["ux"];
    data.uy = var_map["uy"];

    data.number_of_signals = Target_Signals.cols();
    data.number_of_equations = var_map["number_of_equations"];
    data.order_of_equations = var_map["order_of_equations"];


    data.input_connectivity_percentage = var_map["input_connectivity_percentage"];
    data.feedback_connectivity_percentage = 0;

    //data.w_in_initial = -1;
    data.min_input_weight = var_map["min_input_weight"];
    data.max_input_weight = var_map["max_input_weight"];


    data.min_x_position = var_map["min_x_position"];
    data.max_x_position = var_map["max_x_position"];
    data.min_y_position = var_map["min_y_position"];
    data.max_y_position = var_map["max_y_position"];

    data.min_k3 = var_map["min_k3"];
    data.max_k3 = var_map["max_k3"];
    data.min_d3 = var_map["min_d3"];
    data.max_d3 = var_map["max_d3"];

    data.min_k1 = var_map["min_k1"];
    data.max_k1 = var_map["max_k1"];
    data.min_d1 = var_map["min_d1"];
    data.max_d1 = var_map["max_d1"];

    data.dt = var_map["dt"];
    data.t0 = wash_out_time * data.dt;
    data.tmax = (wash_out_time + learning_time + learning_time_test) * data.dt;


    double Mean_Sq = 1000;
    double Mean_Sq_two = 1000;
    double Mean_Sq_three = 1000;

    double MSE = 0;
    double MSE_two = 0;
    double MSE_three = 0;

    double total_MSE = 0;
    double total_MSE_two = 0;
    double total_MSE_three = 0;


    std::vector<double> function_output;
    int number_of_simulations = 1;
    std::mutex m;

    Camera camera(glm::vec3(5.0f, 0.0f, 10.0f));
    bool valid_output = false;

    auto lambda_simulation = [&]() {
        for (int i = 0; i < number_of_simulations; i++) {
            Simulation sim(data, Input_Signals, Target_Signals, wash_out_time, learning_time, learning_time_test, camera, false);
            {
                std::scoped_lock<std::mutex> lock(m);
                /*std::optional<std::vector<double>>opt_learning_matrix = sim.output_LearningMatrix_and_MeanSquaredError();
                if (opt_learning_matrix)
                {*/
                try {
                    function_output = std::move(sim.output_LearningMatrix_and_MeanSquaredError().value());
                }
                catch (const char* msg) {
                    std::cerr << msg << std::endl;
                    return 0;
                }
                //  }

                MSE = std::move(function_output[0]);
                MSE_two = std::move(function_output[1]);
                MSE_three = std::move(function_output[2]);

                total_MSE += MSE;
                total_MSE_two += MSE_two;
                total_MSE_three += MSE_three;

                if (Mean_Sq > MSE && Mean_Sq_two > MSE_two && Mean_Sq_three > MSE_three)
                {
                    Mean_Sq = MSE;
                    std::cout << "The best MSE at the moment is: " << Mean_Sq << "\n";
                    Mean_Sq_two = MSE_two;
                    std::cout << "The best MSE at the moment is: " << Mean_Sq_two << "\n";
                    Mean_Sq_three = MSE_three;
                    std::cout << "The best MSE at the moment is: " << Mean_Sq_three << "\n";
                    sim.output_Output_Signal();

                    valid_output = true;
                }


            }
        }
    };


    std::thread test_thread(lambda_simulation);
    test_thread.join();

    total_MSE = total_MSE / number_of_simulations;
    std::cout << "The average MSE for Volterra: " << total_MSE << "\n";
    std::cout << "The best MSE for Volterra : " << Mean_Sq << "\n";

    total_MSE_two = total_MSE_two / number_of_simulations;
    std::cout << "The average MSE for NARMA: " << total_MSE_two << "\n";
    std::cout << "The best MSE for NARMA : " << Mean_Sq_two << "\n";

    total_MSE_three = total_MSE_three/ number_of_simulations;
    std::cout << "The average MSE for 2ndOrder: " << total_MSE_three << "\n";
    std::cout << "The best MSE for 2ndOrder : " << Mean_Sq_three << "\n";


  
    if (!isnan(total_MSE) && !isnan(total_MSE_two) && !isnan(total_MSE_three) && valid_output)
    {
        ////Plotting Routine goes here
        ////Thanks to pnumerics youtube https://www.youtube.com/channel/UCR6RmA40b4dg4oQw_NoSx6Q
        Gnuplot gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"");

        std::vector<double>  target;
        std::vector<double>  output;

        std::vector<double>  targetTwo;
        std::vector<double>  outputTwo;

        std::vector<double>  targetThree;
        std::vector<double>  outputThree;

        std::vector<double>  outputMerged;

        FileReader target_signal("src/Output/targetsignal.csv", target);
        target_signal.file_read(0);
        FileReader output_signal("src/Output/outputsignal.csv", output);
        output_signal.file_read(0);

        FileReader target_signal_two("src/Output/targetsignal_two.csv", targetTwo);
        target_signal_two.file_read(0);
        FileReader output_signal_two("src/Output/outputsignal_two.csv", outputTwo);
        output_signal_two.file_read(0);

        FileReader target_signal_three("src/Output/targetsignal_three.csv", targetThree);
        target_signal_three.file_read(0);
        FileReader output_signal_three("src/Output/outputsignal_three.csv", outputThree);
        output_signal_three.file_read(0);

        FileReader merged_signal("src/Output/merged_output.csv", outputMerged);
        merged_signal.file_read(0);


        gp << "set multiplot layout 4,1\n";
        gp << "set xlabel 'time [s]'\n";
        gp << "set ylabel 'amplitude'\n";
        gp << "set autoscale\n";
        gp << "set grid\n";

        gp << "set title 'Volterra'\n";
        gp << "unset key\n";
        gp << "plot '-' with lines linestyle 1,"
            << "'-' with lines linestyle 3 title 'VolterraOutput'\n";

        gp.send(target);
        gp.send(output);

        gp << "set title '2ndOrder'\n";
        gp << "unset key\n";
        gp << "plot '-' with lines linestyle 1,"
            << "'-' with lines linestyle 2 title 'NARMAOutput'\n";

        gp.send(targetTwo);
        gp.send(outputTwo);


        gp << "set title 'NARMA'\n";
        gp << "unset key\n";
        gp << "plot '-' with lines linestyle 1,"
            << "'-' with lines linestyle 2 title '2ndOrderOutput'\n";

        gp.send(targetThree);
        gp.send(outputThree);

        gp << "set title 'Merged'\n";
        gp << "unset key\n";
        gp << "plot '-' with lines linestyle 4 \n";

        gp.send(outputMerged);

        gp << "unset multiplot\n";

        std::cin.get();
    }
    else
    {
        std::cout << "ERROR: Simulation has either not finished, encountered an error or the MSE is unsatisfactory \n";
    }
    


}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////Quad and Vanderpol////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void feedback_generator()
{

    InitialDataValues data;
    srand(rdtsc());


    std::vector<double>vanderpol;
    std::vector<double>vanderpol_two;

    std::vector<double>quad;
    std::vector<double>quad_two;

    std::vector<double>lokta_volterra;
    std::vector<double>lokta_volterra_two;

    std::vector<double>input;

    std::cout << "-- Start ---------------------------------------- \n";


    auto read_lambda = [](std::string file_name, std::vector<double>& vec_signal, int column = 1) {
        FileReader input_signal(file_name, vec_signal);
        input_signal.file_read(column);
    };

    std::thread input_thread(read_lambda, "src/Data/inputsignal.csv", std::ref(input));
    input_thread.join();

    std::thread vanderpol_thread(read_lambda, "src/Data/vanderpol.csv", std::ref(vanderpol),0);
    std::thread quad_thread(read_lambda, "src/Data/quad.csv", std::ref(quad),0);
    std::thread lokta_thread(read_lambda,"src/Data/lokta_volterra.csv",std::ref(lokta_volterra),0);
    vanderpol_thread.join();
    quad_thread.join();
    lokta_thread.join();
    std::thread vanderpol_thread_2(read_lambda, "src/Data/vanderpol.csv", std::ref(vanderpol_two),1);
    std::thread quad_thread_2(read_lambda, "src/Data/quad.csv", std::ref(quad_two),1);
    std::thread lokta_thread_2(read_lambda, "src/Data/lokta_volterra.csv", std::ref(lokta_volterra_two), 1);
    vanderpol_thread_2.join();
    quad_thread_2.join();
    lokta_thread_2.join();


    MatrixXd Target_Signals(vanderpol.size(), 6);
    MatrixXd Input_Signals(input.size(), 1);


    Target_Signals.col(0) = Eigen::Map<Eigen::VectorXd>(vanderpol.data(), vanderpol.size());
    Target_Signals.col(1) = Eigen::Map<Eigen::VectorXd>(vanderpol_two.data(), vanderpol_two.size());

    Target_Signals.col(2) = Eigen::Map<Eigen::VectorXd>(quad.data(), quad.size());
    Target_Signals.col(3) = Eigen::Map<Eigen::VectorXd>(quad_two.data(), quad_two.size());

    Target_Signals.col(4) = Eigen::Map<Eigen::VectorXd>(lokta_volterra.data(), lokta_volterra.size());
    Target_Signals.col(5) = Eigen::Map<Eigen::VectorXd>(lokta_volterra_two.data(), lokta_volterra_two.size());

    Input_Signals.col(0) = Eigen::Map<Eigen::VectorXd>(input.data(), input.size());


    std::map <std::string, double>var_map;
    readParameter(var_map, "src/Data/DataFeedback.txt");

    double wash_out_time = var_map["wash_out_time"];
    double learning_time = var_map["learning_time"];
    double learning_time_test = var_map["learning_time_test"];

    // Setting parameters for simulation
    // This should be possible to read in from a text file
    data.N = static_cast<int>(var_map["N"]);
    data.ux = var_map["ux"];
    data.uy = var_map["uy"];

    data.number_of_signals = Target_Signals.cols();
    data.number_of_equations = var_map["number_of_equations"];
    data.order_of_equations = var_map["order_of_equations"];
    

    data.input_connectivity_percentage = var_map["input_connectivity_percentage"];
    data.feedback_connectivity_percentage = var_map["feedback_connectivity_percentage"];

    //data.w_in_initial = -1;
    data.min_input_weight = var_map["min_input_weight"];
    data.max_input_weight = var_map["max_input_weight"];
    data.min_feedback_weight = var_map["min_feedback_weight"];
    data.max_feedback_weight = var_map["max_feedback_weight"];

    data.min_x_position =var_map["min_x_position"];
    data.max_x_position =var_map["max_x_position"];
    data.min_y_position =var_map["min_y_position"];
    data.max_y_position =var_map["max_y_position"];

    data.min_k3 = var_map["min_k3"];
    data.max_k3 = var_map["max_k3"];
    data.min_d3 = var_map["min_d3"];
    data.max_d3 = var_map["max_d3"];

    data.min_k1 = var_map["min_k1"];
    data.max_k1 = var_map["max_k1"];
    data.min_d1 = var_map["min_d1"];
    data.max_d1 = var_map["max_d1"];

    data.dt = var_map["dt"];
    data.t0 = wash_out_time * data.dt;
    data.tmax = (wash_out_time + learning_time + learning_time_test) * data.dt;

    //std::pair<double, double>Mean_Sq;
    double Mean_Sq = 1000;
    double Mean_Sq_two = 1000;
    double Mean_Sq_three = 1000;
    double MSE = 0;
    double MSE_two = 0;
    double MSE_three= 0;
    double total_MSE = 0;
    double total_MSE_two = 0;
    double total_MSE_three = 0;

    std::vector<double>function_output;
    int number_of_simulations = 1;
    std::mutex m;

    Camera camera(glm::vec3(5.0f, 0.0f, 10.0f));
    bool valid_output = false;

    //future make it multihreaded
    auto lambda_simulation = [&]() {

        for (int i = 0; i < number_of_simulations; i++) {
            Simulation sim(data, Input_Signals, Target_Signals, wash_out_time, learning_time, learning_time_test, camera, true);
            {

                std::scoped_lock<std::mutex> lock(m);
                std::vector<double> function_output;
                try {
                    function_output = std::move(sim.output_TestMatrix_and_MeanSquaredError());
                }catch (const char* msg) {
                    std::cerr << msg << std::endl;
                    return 0;
                }
                /*  std::optional<std::vector<double>>opt_learning_matrix =  sim.output_LearningMatrix_and_MeanSquaredError();
                  if (opt_learning_matrix)
                  {*/
                sim.output_LearningMatrix_and_MeanSquaredError();
                //}

                MSE = function_output[0];
                MSE_two = function_output[2];
                MSE_three = function_output[4];

                total_MSE += MSE;
                total_MSE_two += MSE_two;
                total_MSE_three += MSE_three;

                if (Mean_Sq > MSE && Mean_Sq_two > MSE_two && Mean_Sq_three > MSE_three)
                {
                    Mean_Sq = MSE;
                    std::cout << "The best MSE at the moment of Van der Pol is: " << Mean_Sq << "\n";
                    Mean_Sq_two = MSE_two;
                    std::cout << "The best MSE at the moment of Quad is: " << Mean_Sq_two << "\n";
                    Mean_Sq_three = MSE_three;
                    std::cout << "The best MSE at the moment of Lokta-Volterra is: " << Mean_Sq_three << "\n";

                    sim.output_Output_Signal();
                    valid_output = true;
                }
            }
         
        }
    };

    std::thread test_thread(lambda_simulation);
    test_thread.join();

    total_MSE = total_MSE / number_of_simulations;
    std::cout << "The average MSE for Van der Pol: " << total_MSE << "\n";
    std::cout << "The best MSE for Van der Pol : " << Mean_Sq << "\n";

    total_MSE_two = total_MSE_two / number_of_simulations;
    std::cout << "The average MSE for Quad: " << total_MSE_two << "\n";
    std::cout << "The best MSE for Quad : " << Mean_Sq_two << "\n";

    total_MSE_three = total_MSE_three / number_of_simulations;
    std::cout << "The average MSE for Lokta-Volterra: " << total_MSE_three << "\n";
    std::cout << "The best MSE for Lokta-Volterra : " << Mean_Sq_three << "\n";



    if (!isnan(total_MSE) && !isnan(total_MSE_two) && !isnan(total_MSE_three) && valid_output)
    {
        ////Plotting Routine goes here
        ////Thanks to pnumerics youtube https://www.youtube.com/channel/UCR6RmA40b4dg4oQw_NoSx6Q
        Gnuplot gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"");

        std::vector<double>  target;
        std::vector<double>  output;
        std::vector<double>  output_col2;

        std::vector<double>  targetTwo;
        std::vector<double>  outputTwo;
        std::vector<double>  outputTwo_col2;

        std::vector<double>  targetThree;
        std::vector<double>  outputThree;
        std::vector<double>  outputThree_col2;

        std::vector<double>  outputMerged;
        std::vector<double>  targetMerged;
        std::vector<double>  outputMerged_col2;
        std::vector<double>  targetMerged_col2;

        std::vector<double>  feedbackMerged;
        std::vector<double>  outputMerged_Test;

        FileReader target_signal("src/Output/targetsignal.csv", target);
        target_signal.file_read(0);
        FileReader output_signal("src/Output/Results/outputsignal_washout.csv", output);
        output_signal.file_read(0);
        FileReader output_signal_col2("src/Output/Results/outputsignal_washout.csv", output_col2);
        output_signal_col2.file_read(1);

        FileReader target_signal_two("src/Output/targetsignal_two.csv", targetTwo);
        target_signal_two.file_read(0);
        FileReader output_signal_two("src/Output/Results/outputsignal_two_washout.csv", outputTwo);
        output_signal_two.file_read(0);
        FileReader output_signal_two_col2("src/Output/Results/outputsignal_two_washout.csv", outputTwo_col2);
        output_signal_two_col2.file_read(1);

        FileReader target_signal_three("src/Output/targetsignal_three.csv", targetThree);
        target_signal_three.file_read(0);
        FileReader output_signal_three("src/Output/Results/outputsignal_three_washout.csv", outputThree);
        output_signal_three.file_read(0);
        FileReader output_signal_three_col2("src/Output/Results/outputsignal_three_washout.csv", outputThree_col2);
        output_signal_three_col2.file_read(1);

        FileReader merged_signal("src/Output/merged_output.csv", outputMerged);
        merged_signal.file_read(0);
        FileReader merged_target("src/Output/merged_target.csv", targetMerged);
        merged_target.file_read(0);

        FileReader merged_signal_col2("src/Output/merged_output.csv", outputMerged_col2);
        merged_signal_col2.file_read(1);
        FileReader merged_target_col2("src/Output/merged_target.csv", targetMerged_col2);
        merged_target_col2.file_read(1);

        FileReader merged_test_signal("src/Output/merged_test_output.csv", outputMerged_Test);
        merged_test_signal.file_read(0);
        FileReader merged_feedback("src/Output/merged_feedback.csv", feedbackMerged);
        merged_feedback.file_read(0);


        gp << "set multiplot layout 5,1\n";
        gp << "set xlabel 'time [s]'\n";
        gp << "set ylabel 'amplitude'\n";
        gp << "set autoscale\n";
        gp << "set grid\n";

        gp << "set title 'VanderPol'\n";
        gp << "unset key\n";
        gp << "plot '-' with lines linestyle 2,"
            << "'-' with lines linestyle 3 title 'VanderPolOutput'\n";

        gp.send(output);
        gp.send(output_col2);

        gp << "set title 'Quad'\n";
        gp << "unset key\n";
        gp << "plot '-' with lines linestyle 2,"
            << "'-' with lines linestyle 3 title 'QuadOutput'\n";

        gp.send(outputTwo);
        gp.send(outputTwo_col2);

        gp << "set title 'LoktaVolterra'\n";
        gp << "unset key\n";
        gp << "plot '-' with lines linestyle 2,"
            << "'-' with lines linestyle 3 title 'LoktaVolterra'\n";

        gp.send(outputThree);
        gp.send(outputThree_col2);

        gp << "set title 'Merged'\n";
        gp << "unset key\n";
        gp << "plot '-' with lines linestyle 1,"
            << "'-' with lines linestyle 2,"
            << "'-' with lines linestyle 1,"
            << "'-' with lines linestyle 2\n";
        gp.send(outputMerged);
        gp.send(targetMerged);
      
        gp.send(outputMerged_col2);
        gp.send(targetMerged_col2);

        gp << "set title 'MergedTest'\n";
        gp << "unset key\n";
        gp << "plot '-' with lines linestyle 1,"
            << "'-' with lines linestyle 3\n";

        gp.send(feedbackMerged);
        gp.send(outputMerged_Test);

        gp << "unset multiplot\n";

        std::cin.get();
    }
    else
    {
        std::cout << "ERROR: Simulation has either not finished, encountered an error or the MSE is unsatisfactory \n";
    }


   

}


void readParameter(std::map <std::string, double>& var_map,std::string name)
{
    //David Peterson
      //https://stackoverflow.com/questions/27927714/reading-variables-from-a-text-file-c
    std::string line;
    std::string key;
    double value;
    std::ifstream stream(name);
    std::stringstream splitter;

    if (stream)
    {
        while (std::getline(stream, line))
        {
            splitter << line;
            splitter >> key;
            splitter >> value;
            splitter.clear();
            var_map[key] = value;


        }

    }
    else {
        std::cout << "Error Reading File" << std::endl;
    }

}