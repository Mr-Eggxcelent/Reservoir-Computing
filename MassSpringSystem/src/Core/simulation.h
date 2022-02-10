//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Original Code by Alan Quille- Bristol University
//https://github.com/AlanQuille
//The code can be found at https://github.com/AlanQuille/Nonlinear-Mass-Spring-System-BR
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma once
#include <vector>
#include<memory>
#include "system.h"
#include "../Graphics/window.h"

#include <optional>
#include <functional>

#define DEBUG_DRAW 1
#define PERTURBATION 0

extern unsigned int WIDTH;
extern unsigned int HEIGHT;

using namespace Eigen;

class DataSet
{
private:
    double _t0;
    double _dt;
    double _tmax;
    int _maxtimesteps;

public:
    //This loads in the initial values for the signal, t0, tmax and dt
    DataSet(double t0, double tmax, double dt)
    {
        _maxtimesteps = (int)((tmax - t0) / dt);
        _tmax = tmax;
        _dt = dt;
        _t0 = t0;
    }

    //A simple sinewave to test target signal;
    void SineWave(std::vector<double>& sine_wave)
    {
        for (int i = 0; i < _maxtimesteps; i++)
        {
            sine_wave.push_back(sin(_t0 + i * _dt));
            std::cout << sin(_t0 + i * _dt) << "\n";
        }
    }

};


class Simulation
{

private:

    InitialDataValues _data;

    double _t0;
    double _tmax;
    double _dt;

    //These time variables are for the washout, learning phase and test data for the weights calcualated in learning phase.
    unsigned int _wash_out_time;
    unsigned int _learning_time;
    unsigned int _testing_time;

    unsigned int _maxtimesteps;

    bool _feedback_state;

    //This is for the chaos test;
    double _k = 1.000;
    std::string _str = "src/Output/chaoscheck.csv";
    std::string _str2 = "src/Output/LM.csv";


    //This is the learning_weights vector
    std::vector<double> _Learning_Weights;
    //This is the target signal for learning
    MatrixXd _Target_Signal;
    //Input signal for system to be simulated.
    MatrixXd _Input_Signal;

    ////////////////////////////////////////
    //LearningMatrix for learning weight multiplication
    // For collecting data for learning
    //LM is the matrix in total. Lm2 is the matrix for learning and lm3 is the final matrix.

   
    std::vector< std::array< MatrixXd,3>> _LearningMatrix;
    std::vector< std::array< MatrixXd,3>> _Target;
    std::array< MatrixXd, 3>_LearningMatrix_merged;

    //Merged Target
    MatrixXd _mergedTargetSignal;//target signals merged
    MatrixXd _LearningWeightsMatrix;
    MatrixXd _learning_weights;

    MatrixXd _Output;

    std::vector<MatrixXd> _TestMatrix;
    MatrixXd _Merged_Test_Matrix;
    MatrixXd _Test_Feedback;
    MatrixXd _Test_Output;
    std::vector<std::vector<double>> _Test_Feedback_Merged;
    std::vector<std::vector<double>>  _Test_Output_Merged;
    std::vector<std::vector<double>>_Test_Output_Signal;


    //For testing phase.
    VectorXd _Test_Target;
    std::unique_ptr<SpringSystem> _mass_spring;

    //Controls drawing to screen can seperate to a rendering class later
    Window _windowHandle;
    Camera& _camera;

    // timing
    float _deltaTime = 0.0f;
    float _lastFrame = 0.0f;

    SpringSystem::_STATE _state= SpringSystem::_STATE::ORIGINAL;

    //Output

    std::vector<std::vector<double>>_Output_Signal;
    std::vector<std::vector<double>>_Feedback_Signal;
    std::vector<std::vector<double>>_Target_Signal_Vec;

    
    std::vector<std::vector<double>> _OutputMerged;
    std::vector<std::vector<double>> _TargetMerged;
    std::vector<double> _MSE;
    int _number_of_equations;

    bool _key_lock = false;

public:

    //Default constructor
    //Wonder if passing by reference will cause problems while multithreading
    Simulation(InitialDataValues data, MatrixXd Input_Signal, MatrixXd Target_Signals, Camera camera, bool feedback);
    ~Simulation();

    //This is for the csv files to test if there is chaos.
    void input_Magnitude_of_Chaos_Force(double k, const std::string& input, const std::string& input2);

    // Todo: Name is not ideal. Better would be to call it update() or similar
    void execute(bool bias_learning);

    void openLoop();
    void control_feedback();
    void closedLoop();

    void LearningMatrixCreation(unsigned int& index, int index_2, std::array<MatrixXd, 3>& TS, std::array<MatrixXd, 3>& LM, bool draw);
    MatrixXd training_phase(std::array<MatrixXd,3>& LearningMatrix, MatrixXd& y);

    void TestMatrixCreation(unsigned int& index, int index_2, MatrixXd& TM, bool draw);

    //Does the Moore-Penrose pseudoinverse from Eigen library
    void Moore_Penrose_Pseudoinverse(MatrixXd& L);

    //Return Output Signal after learning phase and the mean squared error.
    std::optional<std::vector<double>> output_LearningMatrix_and_MeanSquaredError(bool& success);
    //Test Output
    std::vector<double> output_TestMatrix_and_MeanSquaredError();

    //Outputsignal
    void output_Output_Signal();

    //Return the weights after learning phase as a vector
    std::vector<double>& Return_Learning_Weights();

    //Populate the learning matrix with weights.
    void Populate_Learning_Weights(VectorXd& L);

    //Drawing to screen functions

    void processInput();

    void clearResources();

   
};
