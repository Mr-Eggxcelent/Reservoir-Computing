#pragma once
#include <fstream>
#include<vector>
#include<string>
#include<iostream>

class FileReader
{
public:

	FileReader(std::string file_name,std::vector<double>&vec)
		:_file_name(file_name),_vec(vec)
	{
	}

	void file_read(int column)
	{
		std::ifstream file_input((_file_name).c_str()); file_input.precision(15);
		std::string tmp;

		if (file_input.is_open())
		{
			while (std::getline(file_input, tmp, '\n'))
			{
				_vec_lines.push_back(tmp); //Get each line of the file as a string
			}

			unsigned long s = _vec_lines.size();
			double x;
			for (unsigned int i = 1; i < s; ++i)
			{
				size_t pos = find_column(_vec_lines[i], ",", column); // position of the end of the name of each one in the respective string
				if (pos == 0) {
					x = std::stod(_vec_lines[i].substr(pos, _vec_lines[i].size()));
				}
				else
				{
					x = std::stod(_vec_lines[i].substr(pos + 1, _vec_lines[i].size()));
				}
				x = x * 1;
				_vec.push_back(x); // convert string age to a double
			}
		}
		else
		{
			std::cout << "Error opening file"<<"\n";
		}


		_vec_lines.clear();
	}

	//https://www.oreilly.com/library/view/c-cookbook/0596007612/ch04s11.html
	//Finding the nth Instance of a Substring modified
	size_t find_column(std::string& source, std::string&& pattern, int col)
	{
		if (col == 0) {
			return 0;
		}

		std::string::size_type i = source.find(pattern); // Find the first occurrence

		int j;
		for (j = 1; j < col && i != std::string::npos; ++j)
		{
			i = source.find(pattern, i + 1); // Find the next occurrence
		}
		if (j == col) {
			return(i);
		}
		else {
			return(-1);
		}
	}

private:

	std::string _file_name;
	std::vector<double>& _vec;
	std::vector<std::string> _vec_lines;

};

// Todo: I think it would be useful to define class for reading in files from csv, maybe best even encapsulated in a class
   // there are plenty of these out there
   //  ifstream file_Input ( "/Users/hh14185/Leverhulme_Trust_Proposal_spider_web/Xcode/Nonlinear-Mass-Spring-System-BR/input.csv" );

// TODO:  again best to have that encapsulated in a class - makes the code much cleaner and we will use this code a lot
	//  ifstream file2 ("inputsignal.csv"); // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/

	//  ifstream file_Volterra ("/Users/hh14185/Leverhulme_Trust_Proposal_spider_web/Xcode/Nonlinear-Mass-Spring-System-BR/volterra.csv");