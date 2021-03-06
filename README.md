# Reservoir-Computing
Remember to run it in Release mode.

This project should work for anyone running Visual Studio 2019 on a Windows machine. If not please build the project from the source using the CMake file provided.

The outputs should appear as .csv files in the MassSpringSystem\src\Output\Results folder. There are python files within this folder to plot the output.
Remember to point python to the correct folder when plotting.

If you have not installed GNU plot, make sure to do that. An alternative is to just comment out the internal GNU plotting code in the main file, the program works fine without it.

You can specify the number of simulations you want to carry out within the main.cpp file. If you are happy with the batch of simulations you have run, make sure to go to the Results folder and create a subfolder to store the output. Each time the code is run in this folder the output is overwritten and the old batch will be lost. A specialized save option has yet to be added.

In the simulations.cpp folder there is an option to render the mass-spring system to the screen, switch this on to see the springs in motion. You can't collect data in this mode but it can be helpful if you wish to add some new movement mechanics to the system, and to confirm if their working as expected.

Control the system parameters from the text file within the Data folder.

Whatever C++ file or library that hasn't been written by me has the creators license attached or specified at the top of the file. I'd also like to thank and credit Martin Garrad who I forgot to mention in my thesis for a lot of the code for the mass-spring system. I also didn't mention Daniel Stahlke's library for piping commands from C++ to GNU plot, it was extremely helpful for quickly checking outputs.

Code must be cleaned up, and functions refactored. Offloading calculations to GPU yet to be added.
