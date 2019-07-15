/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                         A.Y. 2015-2016
                  
                Copyright (C) 2017 Stefano Brambilla
======================================================================*/
/*! 
  @file   main.cpp  
  @author Stefano Brambilla <s.brambilla93@gmail.com>   
  @date   September 2016 
  @brief  Main program for test simulations.    
  @details
    We solve the coupled 3D/1D problem of fluid exchange between a 1D 
    network \Lambda and the 3D interstitial tissue \Omega
    
    *****************************************
      Benchmark : single-vessel network 
      Mixed finite elements approximation
      Monolithic resolution by SuperLU 3.0
    *****************************************
    
	See Section "Code verification: test-cases"
 */  
 	 
 	#define M3D1D_VERBOSE_ 
#include <iostream>
#include <problem3d1d.hpp>
#include <transport3d1d.hpp> 
#include <problemHT.hpp>

//! main program
int main(int argc, char *argv[])   
{

	GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
	FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.
 
	try {   
		// Declare a new problem 
		getfem::transport3d1d p; 
		
		// Initialize the problem
                p.problem3d1d::init(argc, argv);
		// Build the monolithic system		
		p.problem3d1d::assembly();
			// Solve the problem
			if(p.problemHT::HEMATOCRIT_TRANSPORT(argc, argv))
				{
				if (!p.problem3d1d::solve()) GMM_ASSERT1(false, "solve procedure has failed");
                                p.problemHT::init(argc, argv);
                                if (!p.problemHT::solve_fixpoint()) GMM_ASSERT1(false, "solve procedure has failed");

                                // Save results in .vtk format
                                p.problemHT::export_vtk();
				
                                /*
                                if(p.OXYGEN_TRANSPORT(argc, argv)) // da aggiungere in transport
                                {
                                    //initialize the transport problem
                                    p.init_transp(argc, argv);
                                    //assemble
                                    p.assembly_transp();
                                    //solve
                                    if (!p.solve_transp()) GMM_ASSERT1(false, "solve procedure has failed");  // the export is in the solve at each time step
                                    // Save results in .vtk format
                                    p.export_vtk();
                                }*/

				}
                        //aggiungere trasporto (con os enza flag)
			else
				{if(!p.problem3d1d::LINEAR_LYMPH())
					{
					// Solve the problem
					if (!p.problemHT::solve_fixpoint()) GMM_ASSERT1(false, "solve procedure has failed");
					
                                        // Add flag with message if oxy-transp YES (and Ht NO)
					}
				else
					{
					// Solve the problem
					if (!p.solve()) GMM_ASSERT1(false, "solve procedure has failed");
					}
				}

		// Save results in .vtk format
                p.problem3d1d::export_vtk();

		// Display some global results: mean pressures, total flow rate
		std::cout << "--- FINAL RESULTS -------------------------" << std::endl; 
		std::cout << "  Pt average            = " << p.problem3d1d::mean_pt()   << std::endl;
		std::cout << "  Pv average            = " << p.problem3d1d::mean_pv()   << std::endl;
		std::cout << "  Network-to-Tissue TFR = " << p.problem3d1d::flow_rate() << std::endl;
		std::cout << "  Lymphatic FR          = " << p.problem3d1d::lymph_flow_rate() << std::endl;
		std::cout << "-------------------------------------------" << std::endl; 	
			      
		   
	}  
      
	GMM_STANDARD_CATCH_ERROR;     
		 
	   
    
		    
		   
	return 0;    
	   
} /* end of main program */   
    
   
  

