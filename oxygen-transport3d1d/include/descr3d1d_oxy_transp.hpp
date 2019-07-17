/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. A.Y. 2018-2019
                  
                Copyright (C) 2018 Riccardo Rosati
======================================================================*/
/*! 
  @file   descr3d1d_oxy_transport3d1d.hpp
  @author Riccardo Rosati <riccardo1.rosati@mail.polimi.it>
  @date   September 2016 - September 2018.
  @brief  Definition of the aux class for algorithm description strings for transport problem.
 */
 

 
/** @defgroup input User-defined parameters  */

#ifndef M3D1D_DESCR3D1D_TRANSP_HPP_
#define M3D1D_DESCR3D1D_TRANSP_HPP_

#include <string>

namespace getfem {

//! Class to import the descriptors of the coupled 3D/1D solver
/*!
	\ingroup input
 */
struct descr3d1d_oxy_transp {

	//General Flags
	//Flag to enable the curve model
	size_type CURVE_PROBLEM;
	//Flag to choose Lymphatic Drainage Curve (0 = sigmoid; 1= linear)
	size_type LINEAR_LYMPHATIC_DRAINAGE;
	//Flag to study the hematocrit distribution all over the network (0 = constant; 1= transport)
	size_type HEMATOCRIT_TRANSPORT;
	//Flag to add advection terms
	size_type ADVECTION;
	//Flag for Type of Viscosity (0 = Vivo - including glycocalyx effect - or 1 = vitro) based on Pries and Secomb works
	size_type Visco_v;
	//Flag to add coupling terms
	size_type COUPLING;
	//Flag to take the stationary problem
	size_type STATIONARY;
	//Flag for the new formulation
	size_type NEW_FORMULATION;
	//Flag to study oxygen transport
	size_type OXYGEN_TRANSPORT;

	//! Absolute path to the vessel mesh file
	std::string MESH_FILEV_OXY;
	//! Identifier of tissue concentration's FEM type
	std::string FEM_TYPET_OT;
	//! Identifier of vessel concentration's FEM type
	std::string FEM_TYPEV_OV;
	//! Identifier of vessel integration method type
	std::string IM_TYPEV_OXY;
	//! Output directory for transport problem
	std::string OUTPUT;


	// Solver information
	//! Identifief of the monolithic solver for transport problem
	std::string SOLVE_METHOD;
	//! Maximum number of iterations (iterative solvers)
	size_type   MAXITER;
	//! Mamimum residual (iterative solvers)
	scalar_type RES; 
	//! Number of target points for the tissue-to-vessel boundary average
	size_type   NInt;
	//! Number of target points for the tissue-to-vessel section average
	size_type   NIntA;

	//! Flag for conforming mesh
	size_type CONFORMING;
	//! Number of 3D region of physical vessel
	size_type SIGMA;
	//! Number of 3D region of tissue outside physical vessel
	size_type OMEGA;
	//! Number of 2D region of physical vessel wall
	size_type GAMMA;
	//! Number of region of the first face of the boundary of the 3d domain
	size_type FACE;
	
	// Utils
	//! File .param
	ftool::md_param FILE_;
	//! Import algorithm specifications from file .param
	void import(ftool::md_param & fname) 
	{
		FILE_ = fname;
		
		MESH_FILEV_OXY  = FILE_.string_value("MESH_FILEV_OXY_TRANSP","1D points file for oxygen transport problem");

		FEM_TYPET_OT   = FILE_.string_value("FEM_TYPET_OT","FEM 3D tissue - concentration");
		FEM_TYPEV_OV   = FILE_.string_value("FEM_TYPEV_OV","FEM 1D vessel - concentration");

		IM_TYPEV_OXY 	= FILE_.string_value("IM_TYPEV_OXY_TRANSP","Name of integration method");

		SOLVE_METHOD_OXY = FILE_.string_value("SOLVE_METHOD_OXY_TRANSP", "Monolithic Solver"); 

		if (SOLVE_METHOD_OXY != "SuperLU") { // iterative solver
			MAXITER  = FILE_.int_value("MAXITER", "Max number of sub-iterations");
			RES = FILE_.real_value("RES"); if (RES == 0.) RES = 2.0e-10;
		}

		NInt = size_type(FILE_.int_value("NInt", "Node numbers on the circle for the nonlocal term"));  
		NIntA = size_type(FILE_.int_value("NIntA", "Node numbers on the radius for the nonlocal term"));  
		OUTPUT = FILE_.string_value("OUTPUT","Output Directory");

		//General Flags
		CURVE_PROBLEM = size_type(FILE_.int_value("CURVE_PROBLEM","Flag to enable the curve model"));
		LINEAR_LYMPHATIC_DRAINAGE = size_type(FILE_.int_value("LINEAR_LYMPHATIC_DRAINAGE","Flag to choose Lymphatic Drainage Curve (0 = sigmoid; 1= linear)"));
		HEMATOCRIT_TRANSPORT = size_type(FILE_.int_value("HEMATOCRIT_TRANSPORT","Flag to study the hematocrit distribution all over the network (0 = constant; 1= transport)"));
		ADVECTION = size_type(FILE_.int_value("ADVECTION","Flag to add advection terms"));
		Visco_v = size_type(FILE_.int_value("Visco_v","Flag for Type of Viscosity (0 = Vivo - including glycocalyx effect - or 1 = vitro) based on Pries and Secomb works"));
		COUPLING = size_type(FILE_.int_value("COUPLING","Flag to add coupling terms"));
		STATIONARY = size_type(FILE_.int_value("STATIONARY","Flag to take the stationary problem"));
		NEW_FORMULATION = size_type(FILE_.int_value("NEW_FORMULATION","Flag for the new formulation"));
		OXYGEN_TRANSPORT = size_type(FILE_.int_value("OXYGEN_TRANSPORT","Flag to study oxygen transport"));

		CONFORMING = size_type(FILE_.int_value("CONFORMING", "Flag for conforming imported mesh by gmsh")); 
		FACE = size_type(FILE_.int_value("FACE", "NNumber of region of the first face of the boundary of the 3d domain"));
		if(CONFORMING)
		{
		SIGMA = size_type(FILE_.int_value("SIGMA", "Number of 3D region of physical vessel"));
		OMEGA = size_type(FILE_.int_value("OMEGA", "Number of 3D region of tissue outside physical vessel"));
		GAMMA = size_type(FILE_.int_value("GAMMA", "Number of 2D region of physical vessel wall"));
		
		}
	}

	//! Overloading of the output operator
	friend std::ostream & operator << (
		std::ostream & out, const descr3d1d_transp & descr
		)
	{ 
		cout << "---- TRANSPORT PROBLEM DESCRIPTORS--------------------------" << endl;
		
		cout << " FEM TYPE  3D concentration     : " << descr.FEM_TYPET_C   << endl;
		cout << " FEM TYPE  1D concentration     : " << descr.FEM_TYPEV_C   << endl;
		cout << "--------------------------------------------------" << endl;

		return out;            
	}

}; /* end of class */

} /* end of namespace */

#endif