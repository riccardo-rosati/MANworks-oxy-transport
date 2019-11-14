/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                         A.Y. 2015-2016
                  
                Copyright (C) 2017 Riccardo Rosati
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
 
#include <iostream>
//#include <AMG_Interface.hpp> //per il MOX ,non mi prendeva i file AMG
#include <problem3d1d.hpp>
#include <oxygen_transport3d1d.hpp> 
#include <problemHT.hpp>

using namespace getfem;

//! main program
int main(int argc, char *argv[])   
{

	GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
	FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.
 
	try {   
		// Declare a new problem 
		oxygen_transport3d1d p; 
		
		// Initialize the problem
                p.problem3d1d::init(argc, argv);
		cout<<"***** Iniziliazzato il problema fluido! *****"<<endl;
		// Build the monolithic system		
		p.problem3d1d::assembly();
		cout<<"***** Assemblato problema fluido! *****"<<endl;
			// Solve the problem
			if(p.problemHT::HEMATOCRIT_TRANSPORT(argc, argv))
				{
				  cout<<"!!C'è trasporto di Ht!!"<<endl;
				if (!p.problem3d1d::solve()) GMM_ASSERT1(false, "solve procedure has failed");
                                p.problemHT::init(argc, argv);
				cout<<"***** Inizializzato l'Ht *****"<<endl;
                                if (!p.problemHT::solve_fixpoint()) GMM_ASSERT1(false, "solve procedure has failed");

        // Display some global results: mean pressures, total flow rate
		std::cout << "--- FINAL RESULTS -------------------------" << std::endl; 
		std::cout << "  Pt average            = " << p.problem3d1d::mean_pt()   << std::endl;
		std::cout << "  Pv average            = " << p.problem3d1d::mean_pv()   << std::endl;
		std::cout << "  Network-to-Tissue TFR = " << p.problem3d1d::flow_rate() << std::endl;
		std::cout << "  Lymphatic FR          = " << p.problem3d1d::lymph_flow_rate() << std::endl;
		std::cout << "-------------------------------------------" << std::endl; 	
		cout<<"***** Risolto il trasporto di Ht *****"<<endl;

                                // Save results in .vtk format
                                p.problemHT::export_vtk();
				
                                
                if(p.OXYGEN_TRANSPORT(argc, argv))
                {
				  cout<<"!!C'è trasporto di ossigeno!!"<<endl;
                                    //initialize the transport problem
                                    p.init_oxy_transp(argc, argv);
				std::cout<<"***** Inizializzato il trasporto di ossigeno *****"<<std::endl;
				
                                    //assemble
                                    p.assembly_oxy_transp();
				std::cout<<"***** Assemblato il trasporto di ossigeno *****"<<std::endl;
                                    //solve
                                    if (!p.solve_oxygen_fixpoint()) GMM_ASSERT1(false, "solve procedure has failed");  // the export is in the solve at each time step
									//if (!p.solve_oxy_transp()) GMM_ASSERT1(false, "solve procedure has failed");  // the export is in the solve at each time step
					cout<<"***** Risolto il trasporto di ossigeno *****"<<endl;
					cout<<"***** Risolvo il trasporto di massa *****"<<endl;
							//p.mass_balance();
                                    // Save results in .vtk format
                                    p.export_vtk_oxy_transp();

				}
				}
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
				p.mass_balance();
		// Save results in .vtk format
                p.problem3d1d::export_vtk();		      
		   
	}  
      
	GMM_STANDARD_CATCH_ERROR;     
		 
	   
    
		    
		  return 0;    
	   
} /* end of main program */   
    
   
  

