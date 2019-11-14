/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. A.Y. 2018-2019
                  
                Copyright (C) 2018 Riccardo Rosati
======================================================================*/
/*! 
  @file   oxygen_transport3d1d.hpp
  @author Riccardo Rosati <riccardo1.rosati@mail.polimi.it>
  @date   September 2016 - September 2018.
  @brief  Declaration of the main class for the 3D/1D coupled transport problem.
 */
 
#ifndef M3D1D_TRANSPORT3D1D_HPP_
#define M3D1D_TRANSPORT3D1D_HPP_


// GetFem++ libraries

#include <gmm/gmm_matrix.h>
#include <getfem/dal_bit_vector.h>

// Project headers
#include <problemHT.hpp>

#include <assembling1d_oxy_transp.hpp>          
#include <assembling3d_oxy_transp.hpp> 
#include <assembling3d1d_oxy_transp.hpp>        
#include <dof3d1d_oxy_transp.hpp>
#include <descr3d1d_oxy_transp.hpp>
#include <param3d1d_oxy_transp.hpp>
#include <utilities_transp.hpp>
#include <node_transp.hpp>
#include "../utilities/muparser/include/muParser.h"


 namespace getfem {

//!	Main class defining the coupled 3D/1D transport problem.
class oxygen_transport3d1d: public problemHT { 

public:
	oxygen_transport3d1d(void) : 
		mf_oxy_Ct(mesht), mf_oxy_Cv(meshv), mf_Ct_Omega(mesht),mf_Ct_Sigma(mesht){} 
	
	// Main methods of class: implement standard and complete transport problem
	//! Initialize the transport problem
	void init_oxy_transp (int argc, char *argv[]);	
	//! Compute the OXYGEN TRANSPORT flag
	bool OXYGEN_TRANSPORT (int argc, char *argv[]);
	//! Assemble the transport problem
	void assembly_oxy_transp (void);
	//! Solve the transport problem
	bool solve_oxy_transp (void);
	//! Solve the FPM for oxygen transport
	bool solve_oxygen_fixpoint (void);
	//! Export the transport solution
	const void export_vtk_oxy_transp (const string & time_suff = "",const string & suff = "");
	//! Compute residuals for mass balance at each junction
	void mass_balance (void);
	//! Compute the adimensional saturation (RR)
	scalar_type dimensioning_saturation (scalar_type cv);
	//! Compute the modified uv in presence of Hemoglobin (RR)
	vector_type modifing_Uvi(vector_type Hi, vector_type uvi, size_type i, vector_type cv);
	//! Compute the oxyhemoglobin concentration before exporting (RR)
	vector_type computing_oxyhemoglobin(vector_type Hi, size_type i, vector_type cv);
	

	//! Getter for solution
	inline vector_type get_UM(void) {return UM_oxy;};

	//Aux methods for interface with problem3d1d class
	//! Initialize the fluid problem
	void init_fluid (int argc, char *argv[]);
	//! Assemble the fluid problem
	void assembly_fluid (void);
	//! Solve the fluid problem
	bool solve_fluid (void);
	//! Export the fluid solution
	const void export_vtk_fluid (const string & suff = "");
	
protected:
	 
	//! Finite Element Method for the tissue oxygen concentration @f$c_t@f$
	mesh_fem mf_oxy_Ct; 
	//! Finite Element Method for the vessel oxygen concentration @f$c_v@f$
	mesh_fem mf_oxy_Cv; 

	//! Finite Element Method for the tissue oxygen concentration @f$c_t@f$
	mesh_fem mf_Ct_Omega; 
	//! Finite Element Method for the tissue oxygen concentration @f$c_t@f$
	mesh_fem mf_Ct_Sigma;

	//! Algorithm description strings (mesh files, FEM types, solver info, ...) 
	descr3d1d_oxy_transp descr_oxy_transp;
	//! Physical parameters
	param3d1d_oxy_transp param_oxy_transp;
	//! Number of degrees of freedom
	dof3d1d_oxy_transp dof_oxy_transp;
		
	
	//! List of BC nodes of the network
	vector< node > BCv_oxy_transp;	
	//! List of BC nodes of the tissue
	vector< node > BCt_oxy_transp;
	//! List of junction nodes of the network
	vector< node_transp > Jv_oxy_transp;

		
	//! Monolithic matrix for the coupled problem
	sparse_matrix_type AM_oxy;
	//! Monolithic array of unknowns for the coupled problem
	vector_type        UM_oxy;
	//! Monolithic right hand side for the coupled problem
	vector_type        FM_oxy;

	// Aux tissue-to-vessel interpolation matrix
	sparse_matrix_type MLIN;
	// Aux tissue-to-vessel average matrix (circumsference)
	sparse_matrix_type MBAR;
	// Aux tissue-to-vessel average matrix (section) 
	sparse_matrix_type MBARBAR;


	// Aux methods for init
	//! Import algorithm specifications
	void import_data_oxy_transp(void);
	//! Import mesh for tissue (3D) and vessel (1D)  
	void build_mesh_oxy_transp(void); 
	//! Set finite elements methods and integration methods 
	void set_im_and_fem_oxy_transp(void);
	//! Build problem parameters
	void build_param_oxy_transp(void);
	//! Build the list of tissue boundary data 
	/*!	Face numbering:
		  0 : {x = 0 }  "back"
		  1 : {x = Lx}  "front"
		  2 : {y = 0 }  "left"
		  3 : {y = Ly}  "right"
		  4 : {z = 0 }  "bottom"
		  5 : {z = Lz}  "top"
	 */
	void build_tissue_boundary_oxy_transp(void);
	//! Build the list of vessel boundary (and junctions) data 
	void build_vessel_boundary_oxy_transp(void);
	
	//Aux method for assembly
	//! Build the monolithic matrix AM_transp by blocks
	void assembly_mat_oxy_transp(void);
	//! Build the monolithic rhs FM_transp by blocks
	void assembly_rhs_oxy_transp(void);
	
	//void update_transp(void);
	//! Aux function for solver: contains the different solving methods and actually solve the system
	bool solver (	const size_type dof1=0, 
			const size_type dof2=0, 
			const size_type dof3=0,
			const size_type dof4=0);	

	//! Set finite elements and integration methods: mf_t should be defined only on Omega_plus	
	void set_im_and_fem_model(void);



	
}; //end of class trasport3d1d


}  //end of namespace

#endif
