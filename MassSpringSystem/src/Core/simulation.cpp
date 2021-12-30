//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Original Code by Alan Quille- Bristol University
//https://github.com/AlanQuille
//The code can be found at https://github.com/AlanQuille/Nonlinear-Mass-Spring-System-BR
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define _USE_MATH_DEFINES
#include <random>
#include <algorithm>
#include <cmath>
#include <ctime>
#include "simulation.h"
#include "Eigen/QR"
#include "Eigen/SVD"
#include "Eigen/LU"
#include <sstream>
#include <fstream>


#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "../Graphics/shader.h"

#include <iostream>
#include <memory>
#include <iomanip>
#include<thread>

#define DEBUG_DRAW 1
void initGL();

using namespace Eigen;


unsigned int WIDTH = 1000;
unsigned int HEIGHT = 1000;

//Aleksandar Haber
//https://github.com/AleksandarHaber/Save-and-Load-Eigen-Cpp-Matrices-Arrays-to-and-from-CSV-files
void saveData(std::string fileName, MatrixXd&  matrix)
{
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

    std::ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}

Simulation::Simulation(InitialDataValues data, MatrixXd IS, MatrixXd TS, int wash_out_time,
    int learning_time, int learning_time_test, Camera camera, bool feedback)
    : _windowHandle(camera, WIDTH, HEIGHT), _camera(camera), _Target_Signal(TS), _Input_Signal(IS)
    , _mass_spring(std::make_unique<SpringSystem>(data, camera, WIDTH, HEIGHT,feedback)), _feedback_state(feedback),_data(data)
{

    _t0 = data.t0;
    _tmax = data.tmax;
    _dt = data.dt;

    _wash_out_time = wash_out_time;
    _learning_time = learning_time;
    _learning_time_test = learning_time_test;

    //Total time
    _maxtimesteps = wash_out_time + learning_time + learning_time_test;
    _number_of_equations = _Target_Signal.cols() / data.order_of_equations;

    
    for (int i = 0; i < _number_of_equations; i++) {

        _LearningMatrix.push_back({ MatrixXd(_maxtimesteps, _mass_spring->get_spring_vec().size()),
                                                  MatrixXd(_learning_time, _mass_spring->get_spring_vec().size()),
                                                  MatrixXd(_learning_time_test, _mass_spring->get_spring_vec().size()) });

        _Target.push_back({ MatrixXd(_maxtimesteps, data.order_of_equations),
                                          MatrixXd(_learning_time, data.order_of_equations),
                                          MatrixXd(_learning_time_test, data.order_of_equations) });

        _TestMatrix.push_back(MatrixXd(_learning_time_test, _mass_spring->get_spring_vec().size()));

    }
    

    _LearningMatrix_merged = { MatrixXd(_maxtimesteps * _number_of_equations,  _LearningMatrix[0][0].cols()),
                               MatrixXd(_learning_time * _number_of_equations,  _LearningMatrix[0][1].cols()),
                               MatrixXd(_learning_time_test * _number_of_equations, _LearningMatrix[0][2].cols()) };


    _mergedTargetSignal = MatrixXd(_Target[0][1].rows() * _number_of_equations, data.order_of_equations);

    _Merged_Test_Matrix = MatrixXd(_learning_time_test * _number_of_equations, _mass_spring->get_spring_vec().size());


    _Output_Signal.resize(_Target_Signal.cols());
    _Target_Signal_Vec.resize(_Target_Signal.cols());
    _Test_Output_Signal.resize(_Target_Signal.cols());
    _Feedback_Signal.resize(_Target_Signal.cols());

    _OutputMerged.resize(data.order_of_equations);
    _TargetMerged.resize(data.order_of_equations);

    _Test_Output_Merged.resize(data.order_of_equations);
    _Test_Feedback_Merged.resize(data.order_of_equations);

    execute(true);

}



Simulation::~Simulation()
{
}



void Simulation::input_Magnitude_of_Chaos_Force(double k, const std::string& input, const std::string& input2)
{
    _k = k;
    _str = input;
    _str2 = input2;
}



void Simulation::execute(bool bias_learning)
{

    /////////////////////////////////
    //  SIMULATION LOOP
    /////////////////////////////////

    PROFILE_FUNCTION();



#if (DEBUG_DRAW==1)

    _windowHandle.initWindow();
    initGL();
    _mass_spring->init_draw();
    unsigned int i = 0;

   
    _state = SpringSystem::_STATE::ORIGINAL;
    while (i < _maxtimesteps && !glfwWindowShouldClose(_windowHandle.getHandle())) {
        // per-frame time logic
       // --------------------
        float currentFrame = glfwGetTime();
        _deltaTime = currentFrame - _lastFrame;
        _lastFrame = currentFrame;
   
        processInput();
   
        glm::mat4 projection = glm::perspective(glm::radians(_camera.Zoom), ((float)_windowHandle.getWidth()) / ((float)_windowHandle.getHeight()), 0.1f, 100.0f);
        glm::mat4 view = _camera.GetViewMatrix();
   
        glClearColor(0, 0, 0, 0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
        
        LearningMatrixCreation(i, 0, _Target[0], _LearningMatrix[0], true);
   
   
        _mass_spring->render(projection, view);
   
        i++;
   
        glfwSwapBuffers(_windowHandle.getHandle());
        glfwPollEvents();
    }
    
    clearResources();
#else

    if (_feedback_state)
    {
        openLoop(); closedLoop();
    }else
    {
        openLoop();
    }

#endif

}



void Simulation::openLoop()
{
    size_t row_offset = 0;
    size_t row_offset_target = 0;


    for (unsigned int i = 0; i < _number_of_equations; i++)
    {
        _state = static_cast<SpringSystem::_STATE>(i);
        for (unsigned int k = 0; k < _wash_out_time; k++)
        {
            for (unsigned int j = 0; j < _mass_spring->get_spring_vec().size(); j++)
            {
                _mass_spring->calculate_forces(j, _dt);
            }
            for (auto& l : _mass_spring->_n)
            {
             //Change the node position, velocity and acceleration in response.
             l.update(_dt);
             _mass_spring->buckle_system(l, _dt, _state, false);
             l.reset_forces();
            }
        }

        switch (_state) {
        case SpringSystem::_STATE::UP: _mass_spring->Output_For_Plot("up");
            break;
        case SpringSystem::_STATE::DOWN:_mass_spring->Output_For_Plot("down");
            break;
        }

       for (unsigned int j = 0; j < _maxtimesteps; j++)
       {
           LearningMatrixCreation(j, (_data.order_of_equations * i), _Target[i], _LearningMatrix[i], false);
       }
    }

    ////Vertical Stacking by tangy
    ////https://stackoverflow.com/questions/21496157/eigen-how-to-concatenate-matrix-along-a-specific-dimension
    for (size_t i = 0; i < _number_of_equations; ++i) {

        for (size_t j = 0; j < 3; ++j) {
            long cur_rows = _LearningMatrix[j][i].rows();//first index is the equation, second is which learning matrix of that specific equation
            _LearningMatrix_merged[i].middleRows(row_offset, cur_rows) = _LearningMatrix[j][i];
            row_offset += cur_rows;
        }

        row_offset = 0;
        long cur_rows = _Target[i][1].rows();
        _mergedTargetSignal.middleRows(row_offset_target, cur_rows) = _Target[i][1];
        row_offset_target += cur_rows;
    }

    
    //Jacobian singular value decomposition for Moore Penrose pseudoinverse
    //Learning test matrix and learning target at end of signal;
    _learning_weights = training_phase(_LearningMatrix_merged, _mergedTargetSignal);
    _Output = _LearningMatrix_merged[2] * _learning_weights;

}


void Simulation::closedLoop()
{
    size_t  row_offset = 0;
    int stride = 0;

    _Test_Feedback = MatrixXd(_Output.rows() / _number_of_equations, _Target_Signal.cols());

    for (int j = 0; j < _Target_Signal.cols(); j++) {

        if (j >= 2 && j < 4) {
            stride = _Output.rows() / _number_of_equations;
        } else if (j >= 4) {
            stride = (_Output.rows() / _number_of_equations) * 2;
        } else {
            stride = 0;
        }

        _Test_Feedback.block(0, j, _Output.rows() / _number_of_equations, 1) = _Output.block(stride, (j % 2 == 0) ? 0 : 1, _Output.rows() / _number_of_equations, 1);
    }


    for (unsigned int i = 0; i < _number_of_equations; i++)
    {
        _state = static_cast<SpringSystem::_STATE>(i);

        for (unsigned int j = 0; j < _Test_Feedback.rows(); j++)
        {
            TestMatrixCreation(j, (_data.order_of_equations * i), _TestMatrix[i], false);
        }

        long cur_rows = _TestMatrix[i].rows();
        _Merged_Test_Matrix.middleRows(row_offset, cur_rows) = _TestMatrix[i];
        row_offset += cur_rows;

    }

   
    //bias_learning
    _Merged_Test_Matrix.conservativeResize(_Merged_Test_Matrix.rows(), _Merged_Test_Matrix.cols() + 1);
    _Merged_Test_Matrix.col(_LearningMatrix_merged[1].cols() - 1) = VectorXd::Ones(_learning_time_test);
    _Test_Output = _Merged_Test_Matrix * _learning_weights;


}

void Simulation::LearningMatrixCreation(unsigned int& index, int index_2,std::array<MatrixXd,3>& TS, std::array<MatrixXd, 3>& LM, bool draw)
{ 
    TS[0].row(index)= _Target_Signal.block(index, index_2, 1, _data.order_of_equations);
    if (index >= _wash_out_time && index < (_wash_out_time + _learning_time)) TS[1].row(index - _wash_out_time) = _Target_Signal.block(index, index_2, 1, _data.order_of_equations);
    if (index >= (_wash_out_time + _learning_time)) TS[2].row(index - _wash_out_time - _learning_time) = _Target_Signal.block(index, index_2, 1, _data.order_of_equations);

    if (_feedback_state) {
        for (unsigned int j = 0; j < _mass_spring->get_spring_vec().size(); j++)
        {
            _mass_spring->calculate_forces(j, _dt);

            LM[0](index, j) = _mass_spring->get_spring_vec()[j].get_length() + Utility::White_Noise_Generator(_feedback_state);  // Todo: update is not needed for target signal
            if (index >= _wash_out_time && index < (_wash_out_time + _learning_time)) LM[1](index - _wash_out_time, j) = _mass_spring->get_spring_vec()[j].get_length() + Utility::White_Noise_Generator(_feedback_state);
            if (index >= (_wash_out_time + _learning_time)) LM[2](index - _wash_out_time - _learning_time, j) = _mass_spring->get_spring_vec()[j].get_length() + Utility::White_Noise_Generator(_feedback_state);
        }
     }else {
        for (unsigned int j = 0; j < _mass_spring->get_spring_vec().size(); j++)
        {
            _mass_spring->calculate_forces(j, _dt);

            LM[0](index, j) = _mass_spring->get_spring_vec()[j].get_length();  // Todo: update is not needed for target signal
            if (index >= _wash_out_time && index < (_wash_out_time + _learning_time)) LM[1](index - _wash_out_time, j) = _mass_spring->get_spring_vec()[j].get_length();
            if (index >= (_wash_out_time + _learning_time)) LM[2](index - _wash_out_time - _learning_time, j) = _mass_spring->get_spring_vec()[j].get_length();
        }
    }

    _mass_spring->update_reset_system(index,index_2, _Input_Signal, _Target_Signal, _dt, _state, draw);

}

void Simulation::TestMatrixCreation(unsigned int& index, int index_2,MatrixXd& TM, bool draw)
{
    for (unsigned int j = 0; j < _mass_spring->get_spring_vec().size(); j++)
    {
        _mass_spring->calculate_forces(j, _dt);
        TM(index, j) = _mass_spring->get_spring_vec()[j].get_length() ;
    }

    _mass_spring->update_reset_system(index, index_2, _Input_Signal, _Test_Feedback, _dt, _state, draw);
}


MatrixXd Simulation::training_phase(std::array<MatrixXd,3>& LearningMatrix, MatrixXd& y)
{
    //bias_learning
    LearningMatrix[1].conservativeResize(LearningMatrix[1].rows(), LearningMatrix[1].cols() + 1);
    LearningMatrix[1].col(LearningMatrix[1].cols() - 1) = VectorXd::Ones(_learning_time);


    JacobiSVD<MatrixXd> svd(LearningMatrix[1], ComputeThinU | ComputeThinV);
    MatrixXd Cp = svd.matrixV() * (svd.singularValues().asDiagonal()).inverse() * svd.matrixU().transpose();
  
    _LearningWeightsMatrix = Cp * y;

    //bias_learning
    LearningMatrix[2].conservativeResize(LearningMatrix[2].rows(),LearningMatrix[2].cols() + 1);
    LearningMatrix[2].col(LearningMatrix[1].cols() - 1) = VectorXd::Ones(_learning_time_test);
    
    return _LearningWeightsMatrix;
}


void Simulation::Moore_Penrose_Pseudoinverse(MatrixXd& L)
{
    L = L.completeOrthogonalDecomposition().pseudoInverse();
}

std::vector<double>& Simulation::Return_Learning_Weights()
{
    return _Learning_Weights;
}


void Simulation::Populate_Learning_Weights(VectorXd& L)
{
    for (unsigned int j = 0; j < _mass_spring->get_spring_vec().size(); j++)
    {
       
        _Learning_Weights.push_back(L(j));
    }
}

// Also, writing data to hard disk during simualtion slow it down a lot.
// Btw. any graphical output (even to the terminal) slows the process down a lot
// However, you could have every 1000 points and update message to show the use the simualtion is still going
// Btw. it is good to have a functionality to switch off any of these things by the user
std::optional<std::vector<double>> Simulation::output_LearningMatrix_and_MeanSquaredError(bool& success)
{

    double wjej = 0;
    double currenttime = 0;
    double currentvalue = 0;

    unsigned int stride = 0;
    
    if (_Output.size()< (_learning_time_test*_number_of_equations)) {
        throw "Matrix is not filled";
    }
   
    for (unsigned int i = 0; i < _maxtimesteps; i++)
    {
        if (i >= (_wash_out_time + _learning_time))
        {
            for (int j = 0; j < _Output_Signal.size(); j++) {

                _Target_Signal_Vec[j].push_back(_Target_Signal(i, j));

                if (_feedback_state) {
                    if (j >= 2 && j < 4) {
                        stride = _learning_time_test;
                    }
                    else if (j >= 4) {
                        stride = _learning_time_test + _learning_time_test;
                    }
                    else {
                        stride = 0;
                    }

                    _Output_Signal[j].push_back(_Output(i - _wash_out_time - _learning_time + stride, (j % 2 == 0) ? 0 : 1));

                } else {

                    _Output_Signal[j].push_back(_Output(i - _wash_out_time - _learning_time + stride));
                    stride += _learning_time_test;

                }

            }
        }

        stride = 0;
    }

    //Change this  code
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    unsigned int total_size = 0;
    for (int i = 0; i < _number_of_equations; i++)
    {
        total_size += _Output_Signal[i].size();
    }
    
    for (int i = 0; i < _data.order_of_equations; i++)
    {
        _OutputMerged[i].reserve(total_size);
        _TargetMerged[i].reserve(total_size);
    }



    if (_feedback_state == true) {

        for (int i = 0; i < _Output_Signal.size(); i++)
        {
            if (i % 2 == 0) {
                _OutputMerged[0].insert(_OutputMerged[0].end(), _Output_Signal[i].begin(), _Output_Signal[i].end());
                _TargetMerged[0].insert(_TargetMerged[0].end(), _Target_Signal_Vec[i].begin(), _Target_Signal_Vec[i].end());
            }
            else
            {
                _OutputMerged[1].insert(_OutputMerged[1].end(), _Output_Signal[i].begin(), _Output_Signal[i].end());
                _TargetMerged[1].insert(_TargetMerged[1].end(), _Target_Signal_Vec[i].begin(), _Target_Signal_Vec[i].end());

            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////


        for (int i = 0; i < _Output_Signal.size(); i++)
        {
            double temp_MSE =Utility::MSE(_Output_Signal[i], _Target_Signal_Vec[i]);
            if (temp_MSE > 1000)
            {
                success = false;
                break;
            }
        }


        return std::nullopt;

    }else{
              
        for (int i = 0; i < _Output_Signal.size(); i++)
        {
            _OutputMerged[0].insert(_OutputMerged[0].end(), _Output_Signal[i].begin(), _Output_Signal[i].end());
            _TargetMerged[0].insert(_TargetMerged[0].end(), _Target_Signal_Vec[i].begin(), _Target_Signal_Vec[i].end()); 
        }


        for (int i = 0; i < _Output_Signal.size(); i++)
        {
            _MSE.push_back(Utility::MSE(_Output_Signal[i], _Target_Signal_Vec[i]));
        }


       std::cout << "Volterra: The mean squared error of the output signal versus the target signal is: " << _MSE[0] << "\n";
       std::cout << "2nd Order: The mean squared error of the output signal versus the target signal is: " << _MSE[1] << "\n";
       std::cout << "NARMA: The mean squared error of the output signal versus the target signal is: " << _MSE[2] << "\n";
    

        return std::optional<std::vector<double>>(_MSE);

    }



}


std::vector<double> Simulation::output_TestMatrix_and_MeanSquaredError()
{
    unsigned int stride = 0;
    std::vector<std::vector<double>>subVecOutput(_Test_Output_Signal.size());
    std::vector<std::vector<double>>subVecTarget(_Test_Output_Signal.size());
    std::ofstream merged_test_output("src/Output/Results/merged_test_output.csv");  merged_test_output.precision(15);
    std::ofstream feedback_signal("src/Output/Results/merged_feedback.csv");  feedback_signal.precision(15);


    //The Feedback signal vector is not very important we are most interested in the direct comparison between the actual target and the test output
    // Either way we save the feedback for possible use
    for (unsigned int i = 0; i < _Test_Feedback.rows(); i++)
    {
            for (int j = 0; j < _Test_Feedback.cols(); j++) {

                if (j >= 2 && j < 4) {
                    stride = _learning_time_test;
                }else if (j >= 4)
                {
                    stride = _learning_time_test + _learning_time_test;
                }else
                {
                    stride = 0;
                }

                _Feedback_Signal[j].push_back(_Test_Feedback(i, j));
                _Test_Output_Signal[j].push_back(_Test_Output(i + stride, (j % 2 == 0) ? 0 : 1));

            }

        stride = 0;

    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    unsigned int total_size = 0;

    for (int i = 0; i < _number_of_equations; i++)
    {
        total_size += _Test_Output_Signal[i].size();
    }


    for (int i = 0; i < _data.order_of_equations; i++)
    {
        _Test_Output_Merged[i].reserve(total_size);
        _Test_Feedback_Merged[i].reserve(total_size);
    }



    for (int i = 0; i < _Test_Output_Signal.size(); i++)
    {
        if (i % 2 == 0) {
            _Test_Output_Merged[0].insert(_Test_Output_Merged[0].end(), _Test_Output_Signal[i].begin(), _Test_Output_Signal[i].end());
            _Test_Feedback_Merged[0].insert(_Test_Feedback_Merged[0].end(), _Feedback_Signal[i].begin(), _Feedback_Signal[i].end());
        }else{
            _Test_Output_Merged[1].insert(_Test_Output_Merged[1].end(), _Test_Output_Signal[i].begin(), _Test_Output_Signal[i].end());
            _Test_Feedback_Merged[1].insert(_Test_Feedback_Merged[1].end(), _Feedback_Signal[i].begin(), _Feedback_Signal[i].end());

        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int i = 0; i < _Test_Output_Merged[0].size(); i++)
    {
        merged_test_output << _Test_Output_Merged[0][i]<<","<< _Test_Output_Merged[1][i];
        merged_test_output << "\n";

        feedback_signal << _Test_Feedback_Merged[0][i]<<","<< _Test_Feedback_Merged[1][i];
        feedback_signal << "\n";
    }

    auto slice = [](std::vector<double>& arr,
        int X, int Y) {


        // Starting and Ending iterators
        auto start = arr.begin() + X;
        auto end = arr.begin() + Y + 1;

        // To store the sliced vector
        std::vector<double> result(Y - X + 1);

        // Copy vector using copy function()
        copy(start, end, result.begin());

        // Return the final sliced vector
        return result;

    };

    
    //Slicing is done just to ensure that we are not taking the MSE when there still might be some disturbance from the buckling
    for (int i = 0; i < _Output_Signal.size(); i++)
    {
   
        if (_Test_Output_Signal[i].size() < _learning_time_test)
        {
            throw "Slicing has failed";
        }

        subVecOutput[i] = slice(_Test_Output_Signal[i],_learning_time_test-15000, _learning_time_test);
        subVecTarget[i] = slice(_Target_Signal_Vec[i], _learning_time_test-15000, _learning_time_test);

    }


    for (int i = 0; i < _Output_Signal.size(); i++)
    {
        _MSE.push_back(Utility::MSE(subVecOutput[i], subVecTarget[i]));
    }

    std::cout << "Van der Pol Signal 1: The mean squared error of the output signal versus the target signal is: " << _MSE[0] << "\n";
    std::cout << "Van der Pol Signal 2: The mean squared error of the output signal versus the target signal is: " << _MSE[1] << "\n";
    std::cout << "Quad Signal 1: The mean squared error of the output signal versus the target signal is: " << _MSE[2] << "\n";
    std::cout << "Quad Signal 2: The mean squared error of the output signal versus the target signal is: " << _MSE[3] << "\n";
    std::cout << "Lissajous Signal 1: The mean squared error of the output signal versus the target signal is: " << _MSE[4] << "\n";
    std::cout << "Lissajous Signal 2: The mean squared error of the output signal versus the target signal is: " << _MSE[5] << "\n";



    return _MSE;//JUST FOR TESTING
}

void Simulation::output_Output_Signal()
{

    std::ofstream learningweights("src/Output/learningweights.csv"); learningweights.precision(15);

    std::ofstream targetsignal("src/Output/Results/targetsignal.csv");  targetsignal.precision(15);
    std::ofstream targetsignal_two("src/Output/Results/targetsignal_two.csv");  targetsignal_two.precision(15);
    std::ofstream targetsignal_three("src/Output/Results/targetsignal_three.csv");  targetsignal_three.precision(15);

    std::ofstream targetsignal_washout("src/Output/Results/targetsignal_washout.csv");  targetsignal_washout.precision(15);
    std::ofstream targetsignal_two_washout("src/Output/Results/targetsignal_two_washout.csv");  targetsignal_two_washout.precision(15);
    std::ofstream targetsignal_three_washout("src/Output/Results/targetsignal_three_washout.csv");  targetsignal_three_washout.precision(15);

    
    std::ofstream outputsignal("src/Output/Results/outputsignal.csv");  outputsignal.precision(15);
    std::ofstream outputsignal_two("src/Output/Results/outputsignal_two.csv");  outputsignal_two.precision(15);
    std::ofstream outputsignal_three("src/Output/Results/outputsignal_three.csv");  outputsignal_three.precision(15);


    std::ofstream outputsignal_washout("src/Output/Results/outputsignal_washout.csv");  outputsignal_washout.precision(15);
    std::ofstream outputsignal_two_washout("src/Output/Results/outputsignal_two_washout.csv");  outputsignal_two_washout.precision(15);
    std::ofstream outputsignal_three_washout("src/Output/Results/outputsignal_three_washout.csv");  outputsignal_three_washout.precision(15);


    std::ofstream merged_output("src/Output/Results/merged_output.csv");  merged_output.precision(15);
    std::ofstream merged_target("src/Output/Results/merged_target.csv");  merged_target.precision(15);

    std::ofstream chaoscheck(_str);

    double wjej = 0;
    double currenttime = 0;
    double currentvalue = 0;
    double average = 0;
    double std = 0;
    double Mean_squared_error = 0;


    if (_feedback_state == true) {

        for (int i = 0; i < _OutputMerged[0].size(); i++)
        {
            merged_output << _OutputMerged[0][i] << "," << _OutputMerged[1][i];
            merged_output << "\n";

            merged_target << _TargetMerged[0][i] << "," << _TargetMerged[1][i];
            merged_target << "\n";
        }


        for (unsigned int i = 0; i < _learning_time_test; i++)
        {


            outputsignal << _Test_Output_Signal[0][i] << "," << _Test_Output_Signal[1][i];
            outputsignal << "\n";

            outputsignal_two << _Test_Output_Signal[2][i] << "," << _Test_Output_Signal[3][i];
            outputsignal_two << "\n";

            outputsignal_three << _Test_Output_Signal[4][i] << "," << _Test_Output_Signal[5][i];
            outputsignal_three << "\n";

            //just a temp solution until it can take arbitrary value

            targetsignal << _Feedback_Signal[0][i] << "," << _Feedback_Signal[1][i];
            targetsignal << "\n";

            targetsignal_two << _Feedback_Signal[2][i] << "," << _Feedback_Signal[3][i];
            targetsignal_two << "\n";

            targetsignal_three << _Feedback_Signal[4][i] << "," << _Feedback_Signal[5][i];
            targetsignal_three << "\n";

            if (i >= _wash_out_time)
            {
                outputsignal_washout << _Test_Output_Signal[0][i] << "," << _Test_Output_Signal[1][i];
                outputsignal_washout << "\n";

                outputsignal_two_washout << _Test_Output_Signal[2][i] << "," << _Test_Output_Signal[3][i];
                outputsignal_two_washout << "\n";

                outputsignal_three_washout << _Test_Output_Signal[4][i] << "," << _Test_Output_Signal[5][i];
                outputsignal_three_washout << "\n";

                //just a temp solution until it can take arbitrary value

                targetsignal_washout << _Feedback_Signal[0][i] << "," << _Feedback_Signal[1][i];
                targetsignal_washout << "\n";

                targetsignal_two_washout << _Feedback_Signal[2][i] << "," << _Feedback_Signal[3][i];
                targetsignal_two_washout << "\n";

                targetsignal_three_washout << _Feedback_Signal[4][i] << "," << _Feedback_Signal[5][i];
                targetsignal_three_washout << "\n";

            }

        }


    } else{

        for (int i = 0; i < _OutputMerged[0].size(); i++)
        {
            merged_output << _OutputMerged[0][i];
            merged_output << "\n";

            merged_target << _TargetMerged[0][i];
            merged_target << "\n";
        }


        for (int i = 0; i < _maxtimesteps; i++)
        {
            if (i >= (_wash_out_time + _learning_time))
            {

                outputsignal << _Output_Signal[0][i - _wash_out_time - _learning_time];
                outputsignal << "\n";

                outputsignal_two << _Output_Signal[1][i - _wash_out_time - _learning_time];
                outputsignal_two << "\n";

                outputsignal_three << _Output_Signal[2][i - _wash_out_time - _learning_time];
                outputsignal_three << "\n";


                targetsignal << _Target_Signal_Vec[0][i - _wash_out_time - _learning_time];
                targetsignal << "\n";

                targetsignal_two << _Target_Signal_Vec[1][i - _wash_out_time - _learning_time];
                targetsignal_two << "\n";

                targetsignal_three << _Target_Signal_Vec[2][i - _wash_out_time - _learning_time];
                targetsignal_three << "\n";
            }

        }
    }

  
    std::ofstream outFile("src/Output/Results/bestMSEEdge.csv");

    std::ofstream outFile_1("src/Output/Results/bestMSENode_orig.csv");
    std::ofstream outFile_2("src/Output/Results/bestMSENode_up.csv");
    std::ofstream outFile_3("src/Output/Results/bestMSENode_down.csv");

    std::ifstream Edges("src/Output/Edges.csv");

    std::ifstream nodesInfo_orig("src/Output/NodeInfo_original.csv");
    std::ifstream nodesInfo_up("src/Output/NodeInfo_up.csv");
    std::ifstream nodesInfo_down("src/Output/NodeInfo_down.csv");
   

    if (outFile && Edges && nodesInfo_orig)
    {
        outFile << "from" << "," << "to" << std::endl;
        outFile << Edges.rdbuf();
 
        outFile_1 << "Node" << "," << "x_pos" << "," << "y_pos" << "," << "type" << std::endl;
        outFile_1 << nodesInfo_orig.rdbuf();

        outFile_2 << "Node" << "," << "x_pos" << "," << "y_pos" << "," << "type" << std::endl;
        outFile_2 << nodesInfo_up.rdbuf();

        outFile_3 << "Node" << "," << "x_pos" << "," << "y_pos" << "," << "type" << std::endl;
        outFile_3 << nodesInfo_down.rdbuf();
    }
}



void Simulation::processInput()
{

    if (glfwGetKey(_windowHandle.getHandle(), GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(_windowHandle.getHandle(), true);
    if (glfwGetKey(_windowHandle.getHandle(), GLFW_KEY_W) == GLFW_PRESS)
        _camera.ProcessKeyboard(FORWARD, _deltaTime);
    if (glfwGetKey(_windowHandle.getHandle(), GLFW_KEY_S) == GLFW_PRESS)
        _camera.ProcessKeyboard(BACKWARD, _deltaTime);
    if (glfwGetKey(_windowHandle.getHandle(), GLFW_KEY_A) == GLFW_PRESS)
        _camera.ProcessKeyboard(LEFT, _deltaTime);
    if (glfwGetKey(_windowHandle.getHandle(), GLFW_KEY_D) == GLFW_PRESS)
        _camera.ProcessKeyboard(RIGHT, _deltaTime);     

    if (glfwGetKey(_windowHandle.getHandle(), GLFW_KEY_1) == GLFW_PRESS)
        _state=SpringSystem::_STATE::ORIGINAL;
    if (glfwGetKey(_windowHandle.getHandle(), GLFW_KEY_2) == GLFW_PRESS)
        _state=SpringSystem::_STATE::UP;
    if (glfwGetKey(_windowHandle.getHandle(), GLFW_KEY_3) == GLFW_PRESS)
        _state=SpringSystem::_STATE::DOWN;

}

void Simulation::clearResources()
{
    _mass_spring->clearResources();
    glfwTerminate();
}


void initGL()
{
    // glad: load all OpenGL function pointers
      // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << "\n";
    }

    //// configure global opengl state
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);      // 4-byte pixel alignment

    // enable /disable features
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_TEXTURE_2D);

    glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glClearColor(0, 0, 0, 0);                   // background color
    glClearStencil(0);                          // clear stencil buffer
    glClearDepth(1.0f);                         // 0 is near, 1 is far
    glDepthFunc(GL_LEQUAL);


}


