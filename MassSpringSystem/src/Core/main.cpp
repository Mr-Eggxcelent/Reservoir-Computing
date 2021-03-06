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
#include "../Graphics/camera_ortho.h"
#include "concurrency.h"

//graphing
#include "../Graphics/gnuplot-iostream.h"


using namespace Eigen;

void no_feedback_generator();
void feedback_generator();
void readParameter(std::string name, InitialDataValues& data, const int& target_cols);

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

    //Uncomment the one you want to run
    //no_feedback_generator();
    feedback_generator();

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "The time it took for the programme to run in total in milliseconds: ";
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms \n";


}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////Volterra,NARMA,2ndOrder////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void no_feedback_generator()
{

    InitialDataValues data;

    std::vector<double>volterra;
    std::vector<double>second_order;
    std::vector<double>narma;
    std::vector<double>input;
    
    std::cout << "-- Start ------------------------------------------------------ \n";

    auto read_lambda = [](std::string file_name, std::vector<double>& vec_signal,int column = 1) {
        FileReader input_signal(file_name, vec_signal);
        input_signal.file_read(column);
    };

    // will change this once a multi threaded queue is added
    std::thread input_thread(read_lambda, "src/Data/input.csv", std::ref(input),0);
    std::thread volterra_thread(read_lambda, "src/Data/volterra.csv", std::ref(volterra),0);
    std::thread second_order_thread(read_lambda, "src/Data/2ndOrder.csv", std::ref(second_order),0);
    std::thread NARMA_thread(read_lambda, "src/Data/NARMA.csv", std::ref(narma),0);

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


    readParameter("src/Data/DataNoFeedback.txt", data, Target_Signals.cols());


    int number_of_simulations = 10;
    bool valid_output = false;
    bool error_thrown = false;

  
    std::ofstream MSE_Results("src/Output/Results/MSE_Results.csv");  MSE_Results.precision(15);
    Camera camera(0, data.max_x_position+1, 0, data.max_y_position+1, glm::vec3(0.0f, 0.0f, 0.0f));
    std::vector<std::array<double, 6>>MSE_storage(number_of_simulations);
    std::mutex m;
    thread_pool tp;

    double Mean_Sq = 2;
    double Mean_Sq_two = 2;
    double Mean_Sq_three = 2;

    std::atomic<double> MSE = 0;
    std::atomic<double> MSE_two = 0;
    std::atomic<double> MSE_three = 0;

    std::atomic<double> total_MSE = 0;
    std::atomic<double> total_MSE_two = 0;
    std::atomic<double> total_MSE_three = 0;
    
    std::atomic_int ongoing_count(0);
    std::atomic_int completed_count(0);


    //terrible code that works somehow will have to fix later
    auto lambda_simulation = [&]() {

            ongoing_count++;

            Simulation sim(data, Input_Signals, Target_Signals, camera, false);
            {
                std::vector<double> function_output;
                bool success = true;

                //Success does nothing for the non feedback case, since we dont need a closed loop phase
                try {
                    if (auto output_matrix = sim.output_LearningMatrix_and_MeanSquaredError(success))
                    {
                        function_output = std::move(output_matrix.value());
                    }
                }
                catch (const char* msg) {
                    std::cerr << msg << std::endl;
                    error_thrown = true;
                    return 0;
                }
                
                MSE = std::move(function_output[0]);
                MSE_two = std::move(function_output[1]);
                MSE_three = std::move(function_output[2]);

                total_MSE += MSE;
                total_MSE_two += MSE_two;
                total_MSE_three += MSE_three;

                {
                    // Not sure if the mutex is needed here
                    std::scoped_lock<std::mutex> lock(m);
                    for (int j = 0; j < function_output.size(); j++)
                    {
                        MSE_storage[completed_count][j] = function_output[j];
                    }


                    if (((Mean_Sq + Mean_Sq_two + Mean_Sq_three) / 3) > ((MSE + MSE_two + MSE_three) / 3))
                    {
                        Mean_Sq = MSE;
                        std::cout <<"Thread "<<std::this_thread::get_id()<<":"<< "The best MSE at the moment is: " << Mean_Sq << "\n";
                        Mean_Sq_two = MSE_two;
                        std::cout <<"Thread "<< std::this_thread::get_id()<<":"<< "The best MSE at the moment is: " << Mean_Sq_two << "\n";
                        Mean_Sq_three = MSE_three;
                        std::cout << "Thread " << std::this_thread::get_id() << ":" << "The best MSE at the moment is: " << Mean_Sq_three << "\n";
                        sim.output_Output_Signal();

                        valid_output = true;
                    }
                }

                std::cout << function_output[0] <<" "<< std::this_thread::get_id()<< std::endl;
            }

            completed_count++;
            ongoing_count--;
    };

    #if (DEBUG_DRAW==1)
        number_of_simulations = 1;
    #endif

    for (uint64_t i(0); i < number_of_simulations; ++i)
    {
        tp.submit(lambda_simulation);
    }

    while (1)
    {

        std::cout <<"\n"<< "ongoing: " << ongoing_count << ", completed: " << completed_count << " / " << number_of_simulations << std::endl;;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        if (completed_count == number_of_simulations || error_thrown==true)
        {
            break;
        }

    }


    total_MSE = total_MSE / number_of_simulations;
    std::cout << "The average MSE for Volterra: " << total_MSE << "\n";
    std::cout << "The best MSE for Volterra : " << Mean_Sq << "\n";


    total_MSE_two = total_MSE_two / number_of_simulations;
    std::cout << "The average MSE for 2ndOrder: " << total_MSE_two << "\n";
    std::cout << "The best MSE for 2ndOrder : " << Mean_Sq_two << "\n";

    total_MSE_three = total_MSE_three / number_of_simulations;
    std::cout << "The average MSE for NARMA: " << total_MSE_three << "\n";
    std::cout << "The best MSE for NARMA : " << Mean_Sq_three << "\n";

    for (int i = 0; i < number_of_simulations; i++)
    {
        for (int j = 0; j < MSE_storage[0].size(); j++)
        {
            MSE_Results << MSE_storage[i][j];
            (j < MSE_storage[0].size() - 1) ? MSE_Results << "," : MSE_Results << "\n";
        }
    }

    //Comment this out if you dont have gnuplot installed
    if (!isnan(Mean_Sq) && !isnan(Mean_Sq_two) && !isnan(Mean_Sq_three) && valid_output)
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

        FileReader target_signal("src/Output/Results/targetsignal.csv", target);
        target_signal.file_read(0);
        FileReader output_signal("src/Output/Results/outputsignal.csv", output);
        output_signal.file_read(0);

        FileReader target_signal_two("src/Output/Results/targetsignal_two.csv", targetTwo);
        target_signal_two.file_read(0);
        FileReader output_signal_two("src/Output/Results/outputsignal_two.csv", outputTwo);
        output_signal_two.file_read(0);

        FileReader target_signal_three("src/Output/Results/targetsignal_three.csv", targetThree);
        target_signal_three.file_read(0);
        FileReader output_signal_three("src/Output/Results/outputsignal_three.csv", outputThree);
        output_signal_three.file_read(0);

        FileReader merged_signal("src/Output/Results/merged_output.csv", outputMerged);
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
            << "'-' with lines linestyle 2 title '2ndOrder'\n";

        gp.send(targetTwo);
        gp.send(outputTwo);


        gp << "set title 'NARMA'\n";
        gp << "unset key\n";
        gp << "plot '-' with lines linestyle 1,"
            << "'-' with lines linestyle 2 title 'NARMA'\n";

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
/////////////////////////////////////////Quad,Vanderpol and Lissajous//////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void feedback_generator()
{

    InitialDataValues data;

    std::vector<double>signal_1;
    std::vector<double>signal_1_two;

    std::vector<double>signal_2;
    std::vector<double>signal_2_two;

    std::vector<double>signal_3;
    std::vector<double>signal_3_two;

    std::vector<double>input;

    std::cout << "-- Start ---------------------------------------- \n";


    auto read_lambda = [](std::string file_name, std::vector<double>& vec_signal, int column = 1) {
        FileReader input_signal(file_name, vec_signal);
        input_signal.file_read(column);
    };

    // will change this once a multi threaded queue is added
    std::thread input_thread(read_lambda, "src/Data/input.csv", std::ref(input));//uneeded in the feedback case but should be made optional
    input_thread.join();

    std::thread signal_1_thread(read_lambda, "src/Data/vanderpol.csv", std::ref(signal_1),0);
    std::thread signal_2_thread(read_lambda, "src/Data/quad.csv", std::ref(signal_2),0);
    std::thread signal_3_thread(read_lambda,"src/Data/lokta_volterra.csv",std::ref(signal_3),0);
    signal_1_thread.join();
    signal_2_thread.join();
    signal_3_thread.join();
    std::thread signal_1_thread_2(read_lambda, "src/Data/vanderpol.csv", std::ref(signal_1_two),1);
    std::thread signal_2_thread_2(read_lambda, "src/Data/quad.csv", std::ref(signal_2_two),1);
    std::thread signal_3_thread_2(read_lambda, "src/Data/lokta_volterra.csv", std::ref(signal_3_two), 1);
    signal_1_thread_2.join();
    signal_2_thread_2.join();
    signal_3_thread_2.join();


    MatrixXd Target_Signals(signal_1.size(), 6);
    MatrixXd Input_Signals(input.size(), 1);


    Target_Signals.col(0) = Eigen::Map<Eigen::VectorXd>(signal_1.data(), signal_1.size());
    Target_Signals.col(1) = Eigen::Map<Eigen::VectorXd>(signal_1_two.data(), signal_1_two.size());

    Target_Signals.col(2) = Eigen::Map<Eigen::VectorXd>(signal_2.data(), signal_2.size());
    Target_Signals.col(3) = Eigen::Map<Eigen::VectorXd>(signal_2_two.data(), signal_2_two.size());

    Target_Signals.col(4) = Eigen::Map<Eigen::VectorXd>(signal_3.data(), signal_3.size());
    Target_Signals.col(5) = Eigen::Map<Eigen::VectorXd>(signal_3_two.data(), signal_3_two.size());

    Input_Signals.col(0) = Eigen::Map<Eigen::VectorXd>(input.data(), input.size());


    readParameter("src/Data/DataFeedback.txt", data, Target_Signals.cols());

    double Mean_Sq = 1;
    double Mean_Sq_two = 1;
    double Mean_Sq_three = 1;


    std::atomic<double> MSE = 0;
    std::atomic<double> MSE_two = 0;
    std::atomic<double> MSE_three= 0;

    std::atomic<double> total_MSE = 0;
    std::atomic<double> total_MSE_two = 0;
    std::atomic<double> total_MSE_three = 0;


    int number_of_simulations = 10;
    bool valid_output = false;
    bool error_thrown = false;

    std::ofstream MSE_Results("src/Output/Results/MSE_Results.csv");  MSE_Results.precision(15);
    Camera camera(0, 10.0f, 0, 10.0f, glm::vec3(0.0f, 0.0f, 0.0f));
    std::mutex m;
    thread_pool tp;

    std::vector<double>function_output;
    std::vector<std::array<double, 6>>MSE_storage(number_of_simulations);

    std::atomic_int ongoing_count(0);
    std::atomic_int completed_count(0);


    //future make it multihreaded
    auto lambda_simulation = [&]() {

        ongoing_count++;

            Simulation sim(data, Input_Signals, Target_Signals, camera, true);
            {

                std::vector<double> function_output;
                bool success = true;


                try {

                    sim.output_LearningMatrix_and_MeanSquaredError(success);
                    function_output = std::move(sim.output_TestMatrix_and_MeanSquaredError());

                }
                catch (const char* msg) {
                    std::cerr << msg << std::endl;
                    error_thrown = true;
                    return 0;
                }

                MSE = (function_output[0] + function_output[1]) / 2;
                MSE_two = (function_output[2] + function_output[3]) / 2;
                MSE_three = (function_output[4] + function_output[5]) / 2;

                total_MSE += MSE;
                total_MSE_two += MSE_two;
                total_MSE_three += MSE_three;

                {
                    std::scoped_lock<std::mutex> lock(m);

                    for (int j = 0; j < function_output.size(); j++)
                    {
                        MSE_storage[completed_count][j] = function_output[j];
                    }

                    if (((Mean_Sq + Mean_Sq_two + Mean_Sq_three) / 3) > ((MSE + MSE_two + MSE_three) / 3) && success)
                    {
                        Mean_Sq = MSE;
                        std::cout << "Thread " << std::this_thread::get_id() << ":" << "The best MSE at the moment is: " << Mean_Sq << "\n";
                        Mean_Sq_two = MSE_two;
                        std::cout << "Thread " << std::this_thread::get_id() << ":" << "The best MSE at the moment is: " << Mean_Sq_two << "\n";
                        Mean_Sq_three = MSE_three;
                        std::cout << "Thread " << std::this_thread::get_id() << ":" << "The best MSE at the moment is: " << Mean_Sq_three << "\n";

                        sim.output_Output_Signal();
                        valid_output = true;
                    }

                }
            }
         
         completed_count++;
         ongoing_count--;
    };


    #if (DEBUG_DRAW==1)
        number_of_simulations = 1;
    #endif

    for (uint64_t i(0); i < number_of_simulations; ++i)
    {
        tp.submit(lambda_simulation);
    }

    while (1)
    {

        std::cout << "\n" << "ongoing: " << ongoing_count << ", completed: " << completed_count << " / " << number_of_simulations << std::endl;;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        if (completed_count == number_of_simulations || error_thrown == true)
        {
            break;
        }

    }

    total_MSE = total_MSE / number_of_simulations;
    std::cout << "The average MSE for Van der Pol: " << total_MSE << "\n";
    std::cout << "The best MSE for Van der Pol : " << Mean_Sq << "\n";

    total_MSE_two = total_MSE_two / number_of_simulations;
    std::cout << "The average MSE for Quad: " << total_MSE_two << "\n";
    std::cout << "The best MSE for Quad : " << Mean_Sq_two << "\n";

    total_MSE_three = total_MSE_three / number_of_simulations;
    std::cout << "The average MSE for Lissajous: " << total_MSE_three << "\n";
    std::cout << "The best MSE for Lissajous : " << Mean_Sq_three << "\n";

    for (int i = 0; i < number_of_simulations; i++)
    {
        for (int j = 0; j < MSE_storage[0].size(); j++)
        {
            MSE_Results << MSE_storage[i][j];
            (j < MSE_storage[0].size() - 1) ? MSE_Results << "," : MSE_Results << "\n";
        }
    }

    //Comment this out if you dont have gnuplot installed
    if (!isnan(Mean_Sq) && !isnan(Mean_Sq_two) && !isnan(Mean_Sq_three) && valid_output)
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

        FileReader target_signal("src/Output/Results/targetsignal.csv", target);
        target_signal.file_read(0);
        FileReader output_signal("src/Output/Results/outputsignal_washout.csv", output);
        output_signal.file_read(0);
        FileReader output_signal_col2("src/Output/Results/outputsignal_washout.csv", output_col2);
        output_signal_col2.file_read(1);

        FileReader target_signal_two("src/Output/Results/targetsignal_two.csv", targetTwo);
        target_signal_two.file_read(0);
        FileReader output_signal_two("src/Output/Results/outputsignal_two_washout.csv", outputTwo);
        output_signal_two.file_read(0);
        FileReader output_signal_two_col2("src/Output/Results/outputsignal_two_washout.csv", outputTwo_col2);
        output_signal_two_col2.file_read(1);

        FileReader target_signal_three("src/Output/Results/targetsignal_three.csv", targetThree);
        target_signal_three.file_read(0);
        FileReader output_signal_three("src/Output/Results/outputsignal_three_washout.csv", outputThree);
        output_signal_three.file_read(0);
        FileReader output_signal_three_col2("src/Output/Results/outputsignal_three_washout.csv", outputThree_col2);
        output_signal_three_col2.file_read(1);

        FileReader merged_signal("src/Output/Results/merged_output.csv", outputMerged);
        merged_signal.file_read(0);
        FileReader merged_target("src/Output/Results/merged_target.csv", targetMerged);
        merged_target.file_read(0);

        FileReader merged_signal_col2("src/Output/Results/merged_output.csv", outputMerged_col2);
        merged_signal_col2.file_read(1);
        FileReader merged_target_col2("src/Output/Results/merged_target.csv", targetMerged_col2);
        merged_target_col2.file_read(1);

        FileReader merged_test_signal("src/Output/Results/merged_test_output.csv", outputMerged_Test);
        merged_test_signal.file_read(0);
        FileReader merged_feedback("src/Output/Results/merged_feedback.csv", feedbackMerged);
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

        gp << "set title 'Lissajous'\n";
        gp << "unset key\n";
        gp << "plot '-' with lines linestyle 2,"
            << "'-' with lines linestyle 3 title 'Lissajous'\n";

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



void readParameter(std::string name, InitialDataValues& data,const int& target_cols)
{
    //David Peterson
      //https://stackoverflow.com/questions/27927714/reading-variables-from-a-text-file-c
    std::string line;
    std::string key;
    double value;
    std::ifstream stream(name);
    std::stringstream splitter;
    std::map <std::string, double>var_map;

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



    data.wash_out_time = var_map["wash_out_time"];
    data.learning_time = var_map["learning_time"];
    data.testing_time = var_map["testing_time"];

    // Setting parameters for simulation
    data.N = static_cast<int>(var_map["N"]);
    data.ux = var_map["ux"];
    data.uy = var_map["uy"];

    data.number_of_signals = target_cols;
    data.number_of_equations = var_map["number_of_equations"];
    data.order_of_equations = var_map["order_of_equations"];


    data.input_connectivity_percentage = var_map["input_connectivity_percentage"];
    data.feedback_connectivity_percentage = var_map["feedback_connectivity_percentage"];

    //data.w_in_initial = -1;
    data.min_input_weight = var_map["min_input_weight"];
    data.max_input_weight = var_map["max_input_weight"];
    data.min_feedback_weight = var_map["min_feedback_weight"];
    data.max_feedback_weight = var_map["max_feedback_weight"];

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
    data.t0 = data.wash_out_time * data.dt;
    data.tmax = (data.wash_out_time + data.learning_time + data.testing_time) * data.dt;
    data.buckling_percentage = var_map["buckling_percentage"];


}