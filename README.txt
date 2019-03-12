
				   Introduction
/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
This is the benchmark used to run the experiments described in the paper "A Comparison of Methods for Gradient Field Estimation on Simplicial Meshes". For both 2D and 3D case, it consists in two executable: one graphic user interface that allows to visualize the error distribution using heat maps, and one batch program to iterate several experiments and encode the results in a .csv file.


				   Requirements
/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
- Cinolib library (https://github.com/mlivesu/cinolib)
- Qt Creator (https://www.qt.io/download)

				   Getting started
/:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

For every executable, you have to type the path to the folder of the Cinolib library in your PC. For example, in order to run the GUI for the 2D case, you have to type such path after "LIBS_PATH" in squared_domain_demo.pro.

For what concerns batch_programs, you can load .obj and .off files from your PC to run the experiments on different meshes from the ones used in the paper. In such case, you will have to type the path to the folder containing such meshes after "DATA_PATH", in both batch_program_squared_domain.pro and batch_program_cubic_domain.pro.

