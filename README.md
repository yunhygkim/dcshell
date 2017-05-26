# Disjoint Convex Shell
#### Yun-hyeong Kim, Zhonghua Xi, Jyh-Ming Lien

Disjoint convex shell(DC-shell) is a set of disjoint convex objects approximating non-convex overapping object sets.
The DC-shell method we proposed provides better approximation than those created by other methods. 
In addition, DC-shell enables faster collision response and realistic fracturing simulation by preventing convex objects from overlapping themselves.

DC-shell implemets the algorithms described in the following paper: 
"Disjoint Convex Shell and its Applications in Mesh Unfolding", SPM 2017, by Yun-hyeong Kim, Zhonghua Xi, and Jyh-Ming Lien. 
([Web Site](http://masc.cs.gmu.edu/wiki/DCShell) / [Paper]() / [Video]() / [BibTex]())

## Description

The provided code constructs disjoint convex objects from overlapping or segmented parts as inputs.
The code can produce DC-shells created by LSF(least-squares fit) heuristic method, SVM and exact volume optimization methods.
To demonstrate the power of DC-shell, we studied how DC-shell can be used in mesh unfolding. 
The nets of polyhedra we used were created by [software tools](http://masc.cs.gmu.edu/wiki/Origami) developed by the [MASC group](http://masc.cs.gmu.edu) at George Mason University. 

## Requirements

This disjoint convex shell program works on Mac OS X. We listed the required program and library.

* Mac OS X 10.9.5 or newer

	The provided code was tested on a MacBook Air, with Mac OSX 10.9.5 and 10.11.6.

* CMake 2.6 or newer

	You can easily install CMake from [Macports](). To install CMake program using Macports, type this in the terminal:

		$ sudo port install cmake

* CGAL library

	The provided code requires only CGAL library. Other libraries are included in the "lib" folder. To install CGAL library using Macports, type this in the terminal:

		$ sudo port install cgal +qt5

## Instructions

* To compile the provided code, please type the commands below in the "cshell" folder, the root directory:

		$ chmod +x gen.sh
		$ ./gen.sh

	These commands compile the code and create a "cshell" file linked to the execution file.

* To run the provided code, please type the command below in the root dirctory:

		$ ./cshell <model_1.obj> <model_2.obj> ... <model_n.obj>

	**NOTE**: Before run the the code, you need to segment a model if you want try to get disjoint convex objects from the model.

	The program takes one or more OBJ files. You can try some examples under the "models" folder. 

* To create figures, tables, and plots  our paper contains, please consult the details in the the [wiki page](https://github.com/yunhkim/dcshell/wiki). 

See detailed description in https://

## Usage
Once you run the cshell program, it would open both a window(left) and a control panel(right). You might see like this:

	$ ./cshell models/yoshi-sep.obj

<img src="./window.jpg" height="400" alt="window"> <img src="./control_panel.jpg" height="400" alt="control_panel">

Here are somethings to try:
1. Press **h** to show the convex hulls.
2. Press **simplify hulls** to simplify the convex hulls.
3. Press **trim hulls** to make the hulls disjoint.
	* Two important parameters are **surface sample density** and **c-svm C**
			
		* Increase **surface sample density** will reduce the number of samples used and speeds up computation
		* **c-svm C** affects both running time and output quality lower. The value increases both computation time and output quality.

4. Once all convex ojects are trimmed and disjoint, then save the resulted objects by clicking **Save Hulls** and then you can use mesh unfolder to unfold each convex hull.

All of the functions in this program are listed in the [wiki page](https://github.com/yunhkim/dcshell/wiki). 

## Models

We 


