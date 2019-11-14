/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. A.Y. 2015-2016
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   param3d1d_transp.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016.
  @brief  Definition of the aux class for physical parameters.
 */
 
#ifndef M3D1D_PARAM3D1D_TRANSP_HPP_
#define M3D1D_PARAM3D1D_TRANSP_HPP_

#include <mesh1d.hpp>    // import_network_radius
#include <utilities.hpp> // compute_radius

namespace getfem {

//! Class to handle the physical parameter of the coupled 3D/1D model
/*!
	\ingroup input
 */
struct param3d1d_oxy_transp {

	// Dimensional physical parameters (microcirc applications)
	//Diffusivity in the tissue [m^2/s]
	scalar_type Dt_;
	//Diffusivity in the vessels [m^2/s]
	scalar_type Dv_;
	//Permeability of the vessel wall [m/s]
	scalar_type Perm_;
	//hydraulic conductivity of the lymphatic wall [s * m^2/kg]
	scalar_type Lp_LF_;
	// surface area of lymphatic vessels per unit volume of tissue [1/m]
	scalar_type SV_;

	
	//Max rate of oxygen metabolization [ml_O2/ml_B/s]
	scalar_type m0_;
	//Tissue Concentration guess [kg/m^3]
	scalar_type Ct_guess_;
	//Vessel Concentration guess [kg/m^3]
	scalar_type Cv_guess_;
	//Partial pressure at half max rate of metabolization [mmHg]
	scalar_type Pm_50_;
	//Solubility of oxygen in the tissue [kg/(m^3*mmHg]
	scalar_type alpha_t_;
	//Hill constant [-]
	scalar_type delta_;
	//Mean Corpuscolar Hemoglobin Concentration [-]
	scalar_type MCHC_;
	//Hufner factor [-];
	scalar_type N_;
	//Oxygen tissue solubility [kg/(m^3*mmHg);
	scalar_type alpha_pl_;
	//Partial Pressure of oxygen at half saturation [mmHg]
	scalar_type Ps_50_;
	//Maximum concentration
	scalar_type C_;


	//Under-relaxation coefficient
	//scalar_type underOXY_;


	// Dimensionless physical parameters (test-cases)
	//Oxygen consumption rate
	scalar_type M0_;
	//Inverse of Peclet number for tissue
	vector_type At_;
	//Inverse of Peclet number for vessel
	vector_type Av_;
	//Damkohler number (metabolism VS diffusion)
	//vector_type Dalpha_;
	//Magnitude of leakage from the capillary bed
	vector_type Y_;
	//lymphatic drainage
	vector_type Q_pl_;	
	

	// Utils
	//! File .param
	ftool::md_param FILE_;
	//! Finite Element Method for tissue data
	getfem::mesh_fem mf_datat_;
	//! Finite Element Method for vessel data
	getfem::mesh_fem mf_datav_;
	// Methods
	//! Build the arrays of dimensionless parameters
	void build(ftool::md_param & fname,
			const getfem::mesh_fem & mf_datat,
			 const getfem::mesh_fem & mf_datav
			) 
	{
		FILE_ = fname;
		mf_datat_ = mf_datat;
		mf_datav_ = mf_datav;
		size_type dof_datat = mf_datat_.nb_dof();
		size_type dof_datav = mf_datav_.nb_dof();
		 
//		bool IMPORT_RADIUS = FILE_.int_value("IMPORT_RADIUS");
		bool NONDIM_PARAM  = FILE_.int_value("TEST_PARAM");
		bool EXPORT_PARAM  = FILE_.int_value("EXPORT_PARAM");
		bool TEST_ANALYTICAL = FILE_.int_value("TEST_ANALYTICAL");
		
		#ifdef M3D1D_VERBOSE_
		cout << "  Assembling dimensionless parameters Dt, Dv, Dalpha, Q_pl ... "   << endl;
		#endif
		if (NONDIM_PARAM) {
			// Import dimensionless params from FILE_
			scalar_type Atval = FILE_.real_value("At"); 
			scalar_type Avval = FILE_.real_value("Av"); 
			//scalar_type Dalphaval = FILE_.real_value("D_alpha"); 
			scalar_type Yval = FILE_.real_value("Y"); 
			scalar_type Q_plval  = FILE_.real_value("Q_pl"); 
			// Fill the data arrays
			 At_.assign(dof_datat,  Atval);
			 Av_.assign(dof_datav,  Avval);
			 //Dalpha_.assign(dof_datat,  Dalphaval);
			 Y_.assign(dof_datav,  Yval);
			 Q_pl_.assign(dof_datat,  Q_plval);
		} 
		else { 
			// Import dimensional params from FILE_
			scalar_type P_  = FILE_.real_value("P", "average interstitial pressure [Pa]"); 
			scalar_type U_  = FILE_.real_value("U", "characteristic flow speed in the capillary bed [m/s]"); 
			scalar_type d_  = FILE_.real_value("d", "characteristic length of the problem [m]"); 
			scalar_type k_  = FILE_.real_value("k", "permeability of the interstitium [m^2]"); 
			scalar_type Lp_ = FILE_.real_value("Lp", "Hydraulic conductivity of the capillary walls [m^2 s/kg]"); 

			if(TEST_ANALYTICAL)
			{
				Dt_   = FILE_.real_value("Dt_test","Test Diffusivity in the tissue [m^2/s]");
			}
			else
			{
				Dt_   = FILE_.real_value("Dt","Diffusivity in the tissue [m^2/s]");
			}

			Dv_   = FILE_.real_value("Dv","Diffusivity in the vessels [m^2/s]");
			
			//Importo i coefficienti dimensionali per l'ossigeno
//RR: MICHEALIS-MENTEN:
			m0_    = FILE_.real_value("m0","Max oxygen consumption rate [ml_O2/ml_tissue/s]");
			Ct_guess_	= FILE_.real_value("Ct_guess","Tissue Concentration guess []");
			Cv_guess_	= FILE_.real_value("Cv_guess","Vessel Concentration guess []");
			Pm_50_	= FILE_.real_value("Pm_50","Partial pressure at half max rate of metabolization [mmHg]");
			alpha_t_ = FILE_.real_value("alpha_t","Solubility of oxygen in the tissue [ml_O2/(ml_Blood*mmHg]");
	
//RR: Parametri OSSIEMOGLOBINA
			MCHC_	= FILE_.real_value("MCHC","Mean Corpuscolar Hematocrit Concentration [-]");
			N_ = FILE_.real_value("N","Hufner factor [-]");
			delta_ = FILE_.real_value("delta","Hill constant [-]");
			Ps_50_ = FILE_.real_value("Ps_50","Partial pressure at half saturation [mmHg]");
			alpha_pl_ = FILE_.real_value("alpha_pl","oxygen solubility in the plasma [kg/(m^3*mmHg)]");

			C_ =  FILE_.real_value("C","maximum concentration [kg/m^3]");
			//underOXY_ = FILE_.real_value("UNDER_RELAXATION_COEFFICIENT_OXY","Under-relaxation coefficient for Transport Solution");


			Perm_ = FILE_.real_value("Perm","Permeability of the capillary walls [m/s]");
			Lp_LF_ = FILE_.real_value("Lp_LF","hydraulic conductivity of the lymphatic wall [s * m^2/kg]");
			SV_ = FILE_.real_value("SV","surface area of lymphatic vessels per unit volume of tissue [1/m]");
			
			// Compute the dimentionless params
			At_.assign(dof_datat, Dt_/d_/U_);
			Av_.assign(dof_datav, Dv_/d_/U_);
			Y_.assign(dof_datav, Perm_/U_);
			Q_pl_.assign(dof_datat,Lp_LF_*SV_*P_*d_/U_);
			M0_ = m0_*d_/U_/C_;				
		}
	


		// Check values
		GMM_ASSERT1(At_[0] != 0, "wrong tissue diffusivity (At>0 required)"); 
		GMM_ASSERT1(Av_[0] != 0, "wrong vessel bed diffusivity (Av>0 required)");
		//if (Q_[0] == 0) cout << "Warning: uncoupled problem (Q=0)" << endl;
		
		if (EXPORT_PARAM){
/*			std::string ODIR = FILE_.string_value("OutputDir","OutputDirectory");
			getfem::vtk_export exp(ODIR+"radius.vtk");
			exp.exporting(mf_datav_);
			exp.write_mesh();
			exp.write_point_data(mf_datav_, R_, "R");
			getfem::vtk_export expQ(ODIR+"conductivity.vtk");
			expQ.exporting(mf_datav_);
			expQ.write_mesh();
			expQ.write_point_data(mf_datav_, Q_, "Q"); */
		}

	}

	//! Get the radius at a given dof
	//inline scalar_type R  (size_type i) { return R_[i];  } const
	//! Get test diffusivity
	inline scalar_type Dv () { return Dv_; } const
	//! Get the tissue diffusivity at a given dof
	inline scalar_type At (size_type i) { return At_[i]; } const
	//! Get the vessel diffusivity at a given dof
	inline scalar_type Av (size_type i) { return Av_[i]; } const
	//! Get the linphatic drainage at a given dof
	inline scalar_type Q_pl  (size_type i) { return Q_pl_[i];  } const
	//! Get the leakage of the capillary bed at a given dof
	inline scalar_type Y  (size_type i) { return Y_[i];  } const
	//! Get the maximum concentration in the vessel to adimensionalise the hemoadvection contribute
	inline scalar_type C () { return C_;  } const
		
	//RR: 
	//! Get the maximum consumption rate
	inline scalar_type M0  () { return M0_;  } const
	//! Get the tissue concentration guess
	inline scalar_type Ct_guess  () { return Ct_guess_;  } const
	//! Get vessel concentration guess
	inline scalar_type Cv_guess  () { return Cv_guess_;  } const
	//! Get the tissue oxygen solubility
	inline scalar_type alpha_t  () { return alpha_t_;  } const
	//! Get the partial pressure at maximum consumption rate
	inline scalar_type Pm_50  () { return Pm_50_;  } const
	//! Get the time mean corpuscola hematocrit concentration
	inline scalar_type MCHC  () { return MCHC_;  } const
	//! Get the Hufner factor
	inline scalar_type N  () { return N_;  } const
	//! Get the plasma oxygen solubility
	inline scalar_type alpha_pl  () { return alpha_pl_;  } const
	//! Get the partial pressure at half saturation
	inline scalar_type Ps_50  () { return Ps_50_;  } const
	//! Get the Hill constant
	inline scalar_type delta  () { return delta_;  } const
	//! Get the under-relaxation coefficient
	//inline scalar_type underOXY  () { return underOXY_; } const
		
	//! Get the radius at a given mesh_region
	//scalar_type R  (const getfem::mesh_im & mim, const size_type rg) 
	//	return compute_radius(mim, mf_datav_, R_, rg);  
	//}
	//! Get the radius
	//vector_type & R (void) { return R_; }
	//! Get the vessel wall permeabilities
	vector_type & Q_pl (void) { return Q_pl_; }
	//! Get the Dahmkholer number 
	//vector_type & Dalpha (void) { return Dalpha_; }
	//! Get the tissue diffusivity 
	vector_type & At (void) { return At_; }
	//! Get the vessel diffusivity 
	vector_type & Av (void) { return Av_; }
	//! Get the leakage of the capillary bed
	vector_type & Y  (void) { return Y_;  } const
	//! Overloading of the output operator
	friend std::ostream & operator << (
		std::ostream & out, const param3d1d_oxy_transp & param
		)
	{ 
		out << "---  DIMENSIONLESS PHYSICAL PARAMS ------" << endl;
		out << "  At     : "                << param.At_[0] << endl; 
		out << "  Av : "                << param.Av_[0] << endl; 
		out << "  Y      : "                << param.Y_[0] << endl; 
		out << "  Q_pl : "                << param.Q_pl_[0] << endl; 
		out << "  m0 : "                << param.M0_ << endl;
		out << "  Pm_50 : "                << param.Pm_50_ << endl;
		out << "  Alpha_t : "                << param.alpha_t_ << endl;
		out << "  Hufner factor : "                << param.N_ << endl;
		out << "  MCHC : "                << param.MCHC_ << endl;
		out << "  Hill Constant : "                << param.delta_ << endl;
		out << "  Ps_50 : "                << param.Ps_50_ << endl;
		out << "  Alpha_pl : "                << param.alpha_pl_ << endl;
		out << "--------------------------" << endl;

		return out;            
	}

}; /* end of class */

} /* end of namespace */

#endif
