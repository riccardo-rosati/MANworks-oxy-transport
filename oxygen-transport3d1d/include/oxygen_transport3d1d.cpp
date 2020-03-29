/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. A.Y. 2018-2019
                  
                Copyright (C) 2018 Riccardo Rosati
======================================================================*/
/*! 
  @file   oxygen_transport3d1d.cpp
  @author Riccardo Rosati <riccardo1.rosati@mail.polimi.it>
  @date   March 2018 - April 2019.
  @brief Definizione della classe principale per risolvere il problema del trasporto di ossigeno
 */
 
 #include <oxygen_transport3d1d.hpp>
 #include <AMG_Interface.hpp>
 #include <cmath>
 #include "gmm/gmm_inoutput.h"
 #include "getfem/getfem_import.h"
 #include <utilities_transp.hpp> //mi permette di usare peclet e vscale e reciprocal


 namespace getfem {


	// Inizializzo il problema
 	void oxygen_transport3d1d::init_oxy_transp (int argc, char *argv[]) 
 	{
 	#ifdef M3D1D_VERBOSE_
 	std::cout << "initialize transport problem..."<<std::endl<<std::endl;
 	#endif

	//PARAM.read_command_line(argc, argv);
	//1. Import data (algorithm specifications, boundary conditions, ...)	
	//import_data_oxy_transp();	
	//2. Import mesh for tissue (3D) and vessel network (1D)
	build_mesh_oxy_transp();	
	//3. Set finite elements and integration methods
	set_im_and_fem_oxy_transp();
 	//4. Build problem parameters
 	build_param_oxy_transp();
	//5. Build the list of tissue boundary data
	build_tissue_boundary_oxy_transp();
	//6. Build the list of tissue boundary (and junction) data
 	build_vessel_boundary_oxy_transp();
 	}; // end of init
	 
	 
	bool oxygen_transport3d1d::OXYGEN_TRANSPORT(int argc, char *argv[])
	{
	PARAM.read_command_line(argc, argv);
	//1. Import data (algorithm specifications, boundary conditions, ...)	
	import_data_oxy_transp();
		 
	return descr_oxy_transp.OXYGEN_TRANSPORT; 
	}; //end of oxygen transport


 	// Aux methods for init
	
	// Import algorithm specifications
	void oxygen_transport3d1d::import_data_oxy_transp(void)
	{
	#ifdef M3D1D_VERBOSE_
	cout << "Importing descriptors for tissue and vessel problems ..." << endl;
	#endif

	descr.import(PARAM);
	descr_oxy_transp.import(PARAM); //da' lo stesso file in input sia per il fluido che per il trasport di oxy?

	#ifdef M3D1D_VERBOSE_
	cout << descr_oxy_transp;
	#endif
	}; //end of import_data_oxy_transp
	
	 
	
	// Import mesh for tissue (3D) and vessel (1D)  
	void oxygen_transport3d1d::build_mesh_oxy_transp(void){

		//In order to have the boundary conditions for the nodes
		//we need to build again the 1D mesh from another pts file
		//! \todo write import_msh_file_transp which take in account both .pts files: modify transport3d1d::init_fluid such that it calls import_msh_file_transp.
		//! \todo branch index should start from 1, not from 0: the junctions will use the notation that +branch is inflow and -branch is outflow. Branch 0 force the user to not give in the .pts file an outflow for first branch (this is usually ok, but loses generality)
		mesht.clear();
			bool test = 0;
		test = PARAM.int_value("TEST_GEOMETRY");
		if(test==0){
			#ifdef M3D1D_VERBOSE_
			cout << "Importing the 3D mesh for the tissue ...  "   << endl;
			#endif
			string st("gmsh:"+descr.MESH_FILET);
			getfem::import_mesh(st,mesht);
				 //import_msh_file(descr.MESH_FILET, mesht);
			}
		else{
			#ifdef M3D1D_VERBOSE_
			cout << "Building the regular 3D mesh for the tissue ...  "   << endl;
			#endif
			string st("GT='" + PARAM.string_value("GT_T") + "'; " +
						   "NSUBDIV=" + PARAM.string_value("NSUBDIV_T") + "; " +  
						   "ORG=" + PARAM.string_value("ORG_T") + "; " +  
						   "SIZES=" + PARAM.string_value("SIZES_T") + "; " +  
						   "NOISED=" + PARAM.string_value("NOISED_T")); 
			#ifdef M3D1D_VERBOSE_		
			cout << "mesht description: " << st << endl;
			#endif
			regular_mesh(problem3d1d::mesht, st);
			//}
	
	#ifdef M3D1D_VERBOSE_
	cout << "Importing the 1D mesh for the vessel (transport problem)... "   << endl;
	#endif
	std::ifstream ifs(descr_oxy_transp.MESH_FILEV_OXY); //legge il file di input di tipo .pts;

	GMM_ASSERT1(ifs.good(), "impossible to read from file " << descr_oxy_transp.MESH_FILEV_OXY);
	meshv.clear();
		
	bool Import=PARAM.int_value("IMPORT_CURVE");
	bool Curve=PARAM.int_value("CURVE_PROBLEM");

	if(Curve && !Import){
		import_pts_file(ifs, meshv, BCv_oxy_transp, nb_vertices, descr.MESH_TYPEV, param);
	}
	
	else if(Import && !Curve){
		GMM_ASSERT1(0,"If you want to import the curvature, you need to enable CURVE_PROBLEM=1");
	}
	else if(Import && Curve){
		std::ifstream ifc(PARAM.string_value("CURVE_FILE","curvature file location"));
		GMM_ASSERT1(ifc.good(), "impossible to read from file " << PARAM.string_value("CURVE_FILE","curvature file location"));
		
		import_pts_file(ifs,ifc, meshv, BCv_oxy_transp, nb_vertices, descr.MESH_TYPEV, param);

		ifc.close();
	} else{
		import_pts_file(ifs, meshv, BCv_oxy_transp, nb_vertices, descr.MESH_TYPEV);
		//da mesh1d.hpp : prende in input il file .pts (ottenuto con std::ifstream ifs(path)), la mesh 1D, la lista delle BC, i numeri di elementi, la tipologia di mesh (GT_PK(1,1))
		//come output: ho numero la lista di BC
	}

	nb_branches = nb_vertices.size();
	ifs.close();
	}
	}; // end of build_mesh_geometry

	// Set finite elements methods and integration methods 
	void oxygen_transport3d1d::set_im_and_fem_oxy_transp(void)
	{

	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEMs for tissue and vessel problems ..." << endl;
	#endif
	
	
	pfem pf_Ot = fem_descriptor(descr_oxy_transp.FEM_TYPET_OT);
	pfem pf_Ov = fem_descriptor(descr_oxy_transp.FEM_TYPEV_OV);

	//RR
	//pfem pf_coeft_oxy = fem_descriptor(descr_oxy_transp.FEM_TYPET_ODATA);

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for tissue ..." << endl;
	#endif
		

	mf_oxy_Ct.set_finite_element(mesht.convex_index(), pf_Ot);
	
	//RR:
	//mf_coeft_oxy.set_finite_element(mesht.convex_index(), pf_coeft_oxy);

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for vessel branches ..." << endl;
	#endif

	mf_oxy_Cv.set_finite_element(meshv.convex_index(), pf_Ov);

	
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEM dimensions for tissue and vessel problems ..." << endl;
	#endif

	dof_oxy_transp.set(mf_oxy_Ct, mf_oxy_Cv);
	#ifdef M3D1D_VERBOSE_
	cout << std::scientific << dof_oxy_transp;
	#endif
	

	// I had to delete the meshes and build them again, because of issues on the boundary conditions:
	// Now i have to build again also the mesh_fem and mesh_im objects.
	mimv.clear();
	mf_Uvi.clear();
	mf_Pv.clear();
	mf_coefv.clear();
	mf_coefvi.clear();

	mimt.clear();
	mf_Ut.clear();
	mf_Pt.clear();
	mf_coeft.clear();



	problem3d1d::set_im_and_fem();

	if(PARAM.int_value("HEMATOCRIT_TRANSPORT")){
	mf_Hi.clear();
	mf_coefh.clear();

	problemHT::set_im_and_fem();
	}

	}; // end of set_im_and_oxy_fem
	
	
	// Build problem parameters
	void oxygen_transport3d1d::build_param_oxy_transp(void)
	{
	
	#ifdef M3D1D_VERBOSE_
	cout << "Building parameters for tissue and vessel problems ..." << endl;
	#endif

	param.build(PARAM, mf_coeft, mf_coefv,mf_coefvi);
	param_oxy_transp.build(PARAM, mf_coeft, mf_coefv);

	#ifdef M3D1D_VERBOSE_
	cout << param_oxy_transp ;
	#endif
	}; // end of build_param_oxy_transp
  
  
  	//Build boundary regions on tissue
	void
	oxygen_transport3d1d::build_tissue_boundary_oxy_transp (void) 
	{
	#ifdef M3D1D_VERBOSE_
	cout << "Building tissue boundary ..." << endl;
	#endif

	size_type face=descr_oxy_transp.FACE;

	BCt_oxy_transp.clear();
	BCt_oxy_transp.reserve(2*DIMT);
	// Parse BC data
	string label_in = PARAM.string_value("BClabel_transp", "Array of tissue boundary labels");
	string value_in = PARAM.string_value("BCvalue_transp", "Array of tissue boundary values");
	vector<string> labels = split(label_in, ' ');
	vector<string> values = split(value_in, ' ');
	GMM_ASSERT1(labels.size()==2*DIMT, "wrong number of BC labels");
	GMM_ASSERT1(values.size()==2*DIMT, "wrong number of BC values");
	for (unsigned f=0; f<2*DIMT; ++f) 
		{
			BCt_oxy_transp.emplace_back(labels[f], std::stof(values[f]), 0, face+f);
			#ifdef M3D1D_VERBOSE_
			cout << "  face " << f << " : " << BCt_oxy_transp.back() << endl;
			#endif
		} 
	

	// Build mesht regions
	size_type xx=0, yy=1, zz=2;

	mesh_region border_faces;
	outer_faces_of_mesh(mesht, border_faces);

	for (mr_visitor i(border_faces); !i.finished(); ++i) 
		{

		assert(i.is_face());

		// Unit outward normal : used to identify faces
		//! \todo Use getfem 5.0's function select_faces_of_normal?
		base_node un = mesht.normal_of_face_of_convex(i.cv(), i.f());
		un /= gmm::vect_norm2(un);

		     if (gmm::abs(un[xx] + 1.0) < 1.0E-7) 	// back
			mesht.region(face+0).add(i.cv(), i.f());
		else if (gmm::abs(un[xx] - 1.0) < 1.0E-7) 	// front
			mesht.region(face+1).add(i.cv(), i.f());
		else if (gmm::abs(un[yy] + 1.0) < 1.0E-7) 	// left
			mesht.region(face+2).add(i.cv(), i.f());
		else if (gmm::abs(un[yy] - 1.0) < 1.0E-7) 	// right
			mesht.region(face+3).add(i.cv(), i.f());
		else if (gmm::abs(un[zz] + 1.0) < 1.0E-7) 	// bottom
			mesht.region(face+4).add(i.cv(), i.f());
		else if (gmm::abs(un[zz] - 1.0) < 1.0E-7) 	// top
			mesht.region(face+5).add(i.cv(), i.f());
				
		}


//PROVA: ciclo sulla regione della mesh #2 (= la faccia left), per vedere quali facce di quali convessi fanno parte della mesh.region(j)
/*
for(size_type j=0; j<6; j++){
		getfem::mesh_region &rg = mesht.region(j);
for (getfem::mr_visitor i(rg); !i.finished(); ++i) {
  cout << "La mesh_region "<<j<<" contains convex " << i.cv();
  if (i.is_face()) cout  << "face " << i.f() << endl;
}
}
*/

//PROVA: lista di tutti i convessi (con indici locali e globali per ogni faccia del convesso)
/*
    dal::bit_vector nn = mesht.convex_index();
    bgeot::size_type i;
    for (i << nn; i != bgeot::size_type(-1); i << nn) {
      cout << "\n Convex of index " << i << endl;
      bgeot::pconvex_structure cvs = mesht.structure_of_convex(i);
      cout << "\n Number of vertices: " << cvs->nb_points() << endl;
      cout << "\n Number of faces: " << cvs->nb_faces() << endl;
      for (bgeot::short_type f = 0; f < cvs->nb_faces(); ++f) {
        cout << "face " << f << " has " << cvs->nb_points_of_face(f);
        cout << " vertices with local indexes: ";
        for (bgeot::size_type k = 0; k < cvs->nb_points_of_face(f); ++k)
          cout << cvs->ind_points_of_face(f)[k] << " ";
        cout << " and global indexes: ";
        for (bgeot::size_type k = 0; k < cvs->nb_points_of_face(f); ++k)
          cout << mesht.ind_points_of_convex(i)[cvs->ind_points_of_face(f)[k]] << " ";
      }
    }
*/		

}; //end of build_tissue_boundary_transp

	//Build boundary regions on network

	void 
	oxygen_transport3d1d::build_vessel_boundary_oxy_transp(void)
	{
	#ifdef M3D1D_VERBOSE_
	cout << "Building vessel boundary ..." << endl;
	#endif
	try {

	dal::bit_vector junctions; // global idx of junctions vertices in meshv
	dal::bit_vector extrema;   // global idx of extreme vertices in meshv

	Jv_oxy_transp.clear();
	nb_extrema=0; 
	nb_junctions=0;
	
	size_type fer = nb_branches; // first empty region
	GMM_ASSERT1(meshv.has_region(fer)==0, 
		"Overload in meshv region assembling!");

	// List all the convexes
	dal::bit_vector nn = meshv.convex_index();

	bgeot::size_type cv;
	for (cv << nn; cv != bgeot::size_type(-1); cv << nn) {

		bgeot::pconvex_structure cvs = meshv.structure_of_convex(cv);
		if (cvs->nb_points()>2) 
			cerr << "Error: convex #" << cv << "has more than 2 vertices!" << endl;
		if (cvs->nb_faces()>2)  
			cerr << "Error: convex #" << cv << "has more than 2 faces!" << endl;

		// Build regions for BCs and junctions
		// Global idx of mesh vertices
		
		size_type i0 = meshv.ind_points_of_convex(cv)[cvs->ind_points_of_face(1)[0]];
		size_type i1 = meshv.ind_points_of_convex(cv)[cvs->ind_points_of_face(0)[0]];

		// Identify vertex type
		if (meshv.convex_to_point(i0).size()==1){ /* inflow extremum */
			// Update information
			extrema.add(i0);
			nb_extrema++;
			// Build a new region made by a single face
			GMM_ASSERT1(meshv.has_region(fer)==0, 
				"Overload in meshv region assembling!");
			meshv.region(fer).add(cv, 1);
			// Store the current index and then update it
			size_type bc = 0; 
			bool found = false;
			while (!found && (bc<BCv_oxy_transp.size())) {
				found = (i0 == BCv_oxy_transp[bc].idx);
				if (!found) bc++;
			}
			GMM_ASSERT1(found=true, "Miss a boundary node in BCv list!");

			BCv_oxy_transp[bc].rg = fer;
			fer++;
			// Store the containing branch index
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {

				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");
			BCv_oxy_transp[bc].branches.emplace_back(branch); 
		}
		else if (meshv.convex_to_point(i0).size()==2){ /* trivial inflow junction */
			// DO NOTHING

		}
		else if (meshv.convex_to_point(i0).size()>=2){ /* non-trivial inflow junction */
		// Check if junction has been already stored, 
			// if not add to the junction list (J) and build a new region

	//! \todo add multiple times the junction node to the junction region. Generic_assembly looks (apparently) only at indexes of the nodes, not at his coordinates; in this way, when I build the region with the junction node from a certain branch, generic_assembly will not recognize the same node from another branch (probably i look at the basis buildt only on the first branch). In order to use the generic_assembly for junction nodes I should add all the basis to the region (e.g. the same node from all the branches)

			dal::bit_vector tmp; tmp.add(i0);
			if(!junctions.contains(tmp)){
				// Store the junction vertex
				junctions.add(i0);
				nb_junctions++;
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				// Build a new region with idx "first empty region"
				meshv.region(fer).add(cv, 1); // single-face region
				// Create a new junction node
				Jv_oxy_transp.emplace_back("JUN", 0, i0, fer);
				fer++;
			}
			// Search for index of containing branch (\mathcal{P}^{in}_j)
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");
			// Add the inflow branch (to the right junction node)
			size_type jj = 0;
			bool found = false;
			while (!found && jj < nb_junctions){
				found = (i0 == Jv_oxy_transp[jj].idx);
				if (!found) jj++;
			}
			//cout << "Branch -" << branch << " added to junction " << jj << endl;
			Jv_oxy_transp[jj].value += param.R(mimv, branch);
			Jv_oxy_transp[jj].branches.emplace_back(-branch);
			GMM_ASSERT1(branch>0, 
				"Error in network labeling: -0 makes no sense");
		}
		

		
		if (meshv.convex_to_point(i1).size()==1){ 
			size_type bc = 0; 
			bool found = false;
			while (!found && (bc<BCv_oxy_transp.size())) {
				found = (i1 == BCv_oxy_transp[bc].idx);
				if (!found) bc++;
			}
			if (found){ /* outlow extremum */
				extrema.add(i1); 
				nb_extrema++; 
				// Build a new region made by a single face
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				meshv.region(fer).add(cv, 0);
				// Store the current index and then update it
				BCv_oxy_transp[bc].value *= +1.0;
				BCv_oxy_transp[bc].rg = fer; 
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv_oxy_transp[bc].branches.emplace_back(branch); 
			}

			else { // interior -> Mixed point 
				// "MIX" label via post-processing
			// Build a new region made by a single face
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				meshv.region(fer).add(cv, 0);
				BCv_oxy_transp.emplace_back("MIX", 0.0, i1, fer);
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv_oxy_transp.back().branches.emplace_back(branch); 
			}

		}
		else if (meshv.convex_to_point(i1).size()==2){ /* trivial outflow junction */

			// Search for index of first containing branch (\mathcal{P}^{out}_j)
			size_type firstbranch = 0; 
			bool contained = false;
			while (!contained && firstbranch<nb_branches ) {
				contained = meshv.region(firstbranch).is_in(cv);
				if (!contained) firstbranch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i1!");

			// Check if i1 is a trivial junction (or a INT point)
			size_type cv1 = meshv.convex_to_point(i1)[0];
			size_type cv2 = meshv.convex_to_point(i1)[1];
			bool is_junc = (meshv.region(firstbranch).is_in(cv1) < 1 ||
							meshv.region(firstbranch).is_in(cv2) < 1 );
							
			if (is_junc){
				cout << "Found a trivial junction at i1 = " << i1 << endl;
				// Check if jucntion has been already stored, 
				// if not add to the junction list (J) and build a new region
				dal::bit_vector tmp; tmp.add(i1);
				if(!junctions.contains(tmp)){
					// Store the junction vertex
					junctions.add(i1);
					nb_junctions++;
					GMM_ASSERT1(meshv.has_region(fer)==0, 
						"Overload in meshv region assembling!");
					// Build a new region with idx "first empty region"
					meshv.region(fer).add(cv, 0);
					// Create a new junction node
					Jv_oxy_transp.emplace_back("JUN", 0, i1, fer);
					fer++;
				// Search for index of second containing branch (\mathcal{P}^{out}_j)
				size_type secondbranch = 0; 
				size_type secondcv = (( cv1 == cv) ? cv2 : cv1);
				size_type firstcv = (( cv1 != cv) ? cv2 : cv1);
				contained = false;
				while (!contained && secondbranch<nb_branches ) {
					if (secondbranch!=firstbranch)
					contained = meshv.region(secondbranch).is_in(secondcv);
					if (!contained) secondbranch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				// Add the two branches
				scalar_type in;
				in=0;
				if (meshv.ind_points_of_convex(firstcv)[0]==i1) in=-1;
				else if (meshv.ind_points_of_convex(firstcv)[1]==i1) in=+1;
				GMM_ASSERT1(in!=0, "There's something wrong in firstbranch convex index");
				Jv_oxy_transp.back().branches.emplace_back(in*firstbranch);

				in=0;
				if (meshv.ind_points_of_convex(secondcv)[0]==i1) in=-1;
				else if (meshv.ind_points_of_convex(secondcv)[1]==i1) in=+1;
				GMM_ASSERT1(in!=0, "There's something wrong in secondbranch convex index");
				Jv_oxy_transp.back().branches.emplace_back(in*secondbranch);
				Jv_oxy_transp.back().value += param.R(mimv, firstbranch);
				Jv_oxy_transp.back().value += param.R(mimv, secondbranch);
				}
			}
		}
		else if (meshv.convex_to_point(i1).size()>=2){ /* non-trivial outflow junction */

			// Search for index of containing branch (\mathcal{P}^{out}_j)
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");

			// Check if jucntion has been already stored, 
			// if not add to the junction list (J) and build a new region
			dal::bit_vector tmp; tmp.add(i1);
			if(!junctions.contains(tmp)){
				// Store the junction vertex
				junctions.add(i1);
				nb_junctions++;
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				// Build a new region with idx "first empty region"
				meshv.region(fer).add(cv, 0);
				// Create a new junction node
				Jv_oxy_transp.emplace_back("JUN", 0, i1, fer);
				// Add the outflow branch
				Jv_oxy_transp.back().branches.emplace_back(+branch);
				Jv_oxy_transp.back().value += param.R(mimv, branch);
				//cout << "Branch " << branch << " added to junction " << i1 << endl;
				fer++;
			}
			else {
				// Add the outflow branch (to the right junction node)
				size_type jj = 0;
				bool found = false;
				while (!found && jj < nb_junctions){
					found = (i1 == Jv_oxy_transp[jj].idx);
					if (!found) jj++;
				}
				Jv_oxy_transp[jj].branches.emplace_back(+branch);
				Jv_oxy_transp[jj].value += param.R(mimv, branch);
				//cout << "Branch " << branch << " added to junction " << jj << endl;
			}
		}

	} /* end of convexes loop */  

	// Ckeck network assembly
	#ifdef M3D1D_VERBOSE_
	cout << "--- NETWORK ASSEMBLY ------------------ "   << endl;
	cout << "  Branches:   " << nb_branches << endl
		 << "  Vertices:   " << nn.size()+1 << endl;
	cout << "  Extrema:    " << extrema << endl;	  
	for (size_type i=0; i<BCv_oxy_transp.size(); ++i)
		cout << "    -  label=" << BCv_oxy_transp[i].label 
			 << ", value=" << BCv_oxy_transp[i].value << ", ind=" << BCv_oxy_transp[i].idx 
			 << ", rg=" << BCv_oxy_transp[i].rg << ", branches=" << BCv_oxy_transp[i].branches << endl; 
	cout << "  Junctions: " << junctions << endl;
	for (size_type i=0; i<Jv_oxy_transp.size(); ++i)
		cout << "    -  label=" << Jv_oxy_transp[i].label 
			 << ", value=" << Jv_oxy_transp[i].value << ", ind=" << Jv_oxy_transp[i].idx 
			 << ", rg=" << Jv_oxy_transp[i].rg << ", branches=" << Jv_oxy_transp[i].branches << endl; 
	cout<< "---------------------------------------- "<<endl;
	#endif

	} 
	GMM_STANDARD_CATCH_ERROR; // catches standard errors

	//NB: Getfem, per ogni branch mette i convessi/elementi in ordine, ma globalmente gli inidici/vertici dei convessi non lo sono: saltano, se uso ind_points_of_convex(mrv.cv())
	//invece vengono scelti correttamente con ind_basic_dof_of_element(mrv.cv()) --> ridaà un aray con gli inidici globali dei vertici dell'elemento mrv.cv()

	} /* end of build_vessel_boundary_oxy_transp */

	//////////////////////////// ASSEMBLAGGIO ////////////////////////////
vector_type 
oxygen_transport3d1d::modifing_Uvi(vector_type Hi, vector_type uvi, size_type i, vector_type Cv){
	/* prende in input il vettore ematocrito e il vettore velocità per il branch i-esimo e la concentrazione per tutti vasi.
		la funzione fa in modo tale da calcolare psi, per ogni branch, e successivamente moltiplica la velocità di quel branch
		per (1+psi) --> uv(1+psi), così da considerare la presenza dell'emoglobina  

	*/
		vector_type new_uvi(mf_Uvi[i].nb_dof()); gmm::clear(new_uvi);
		vector_type psi(mf_Hi[i].nb_dof()); gmm::clear(psi); //salvo i valori di PSI = ossigeno legato all'emoglobina
		vector_type new_psi (mf_Uvi[i].nb_dof()); gmm::clear(new_psi);//la nuova PSI: (1+PSI) dopo interpola

		scalar_type k1;
		k1 = param_oxy_transp.MCHC() * param_oxy_transp.N() / param_oxy_transp.C();

		mesh_region &rg_branch = meshv.region(i);

		//salvo i dof della mesh di Cv e Uv e Ht
		vector_type dof_enum_C, dof_enum_U, dof_enum_H, dof_enum_PSI;

size_type pos=0;
	for(getfem::mr_visitor mrv(rg_branch); !mrv.finished(); ++mrv){
		//cout<<"\n Elemento della mesh fem dell'ematocrito: "<<mrv.cv()<<endl;
		if(pos==0){
			dof_enum_H.emplace_back(mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0]); pos++;
			dof_enum_H.emplace_back(mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]); pos++;
			/*ht_i[pos] = Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0]];
			cout << " ht_i [" << pos << "] = Hi[" << mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0] << "]" <<"="<<Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0]]<< endl;
			pos ++;
			ht_i[pos] = Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " ht_i [" << pos << "] = Hi[" << mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;*/
			}
			else{
			dof_enum_H.emplace_back(mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]); pos ++;
			/*ht_i[pos] = Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " ht_i [" << pos << "] = Hi[" << mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]]<< endl;*/
			}
	}

//cout<<"\n Pos= "<<pos<<"\n Branch #"<<i<<": DOF per la mesh dell'ematocrito= "<<dof_enum_H<<endl;
//	cout<<"***************"<<endl;
	pos=0;
	for(getfem::mr_visitor mrv(rg_branch); !mrv.finished(); ++mrv){
		//cout<<"\n Elemento della mesh fem della velocità: "<<mrv.cv()<<endl;
		if(pos==0){
			dof_enum_U.emplace_back(mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[0]); pos++;
			dof_enum_U.emplace_back(mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[1]); pos++;
			dof_enum_U.emplace_back(mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[2]); pos++;
			/*uv_i[pos] = uvi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[0]];
			cout << " Uvi [" << pos << "] = uvi[" << mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[0] << "]" <<"="<<uvi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[0]]<< endl;
			pos ++;
			uv_i[pos] = uvi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " Uvi [" << pos << "] = uvi[" << mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<uvi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;
			uv_i[pos] = uvi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[2]];
			cout << " Uvi [" << pos << "] = uvi[" << mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[2] << "]" <<"="<<uvi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[2]]<< endl;
			pos ++;*/
			}
		else{
			dof_enum_U.emplace_back(mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[1]); pos++;
			dof_enum_U.emplace_back(mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[2]); pos++;
			/*uv_i[pos] = uvi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " Uvi [" << pos << "] = uvi[" << mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<uvi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;
			uv_i[pos] = uvi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[2]];
			cout << " Uvi [" << pos << "] = uvi[" << mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[2] << "]" <<"="<<uvi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[2]]<< endl;
			pos ++;*/
			}
	}
//cout<<"\n Pos= "<<pos<<"\n Branch #"<<i<<": DOF per la mesh della veclotià= "<<dof_enum_U<<endl;
//	cout<<"***********************"<<endl;
		pos=0;
	for (getfem::mr_visitor mrv(rg_branch); !mrv.finished(); ++mrv)
		{
		//cout<<"\n Elemento della mesh fem della concentrazione: "<<mrv.cv()<<endl;
		if(pos == 0){
			dof_enum_C.emplace_back(mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0]); pos++;
			dof_enum_C.emplace_back(mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]); pos++;
			/*cv_i[pos] = Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0]];
			cout << " cv_i [" << pos << "] = cv[" << mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0] << "]" <<"="<<Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0]]<< endl;
			pos ++;
			cv_i[pos] = Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " cv_i [" << pos << "] = cv[" << mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;*/
			}
		else{
			dof_enum_C.emplace_back(mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]); pos++;;
			/*cv_i[pos] = Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " cv_i [" << pos << "] = cv[" << mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;*/
			}
		}
//cout<<"\n Pos= "<<pos<<"\n Branch #"<<i<<": DOF per la mesh della concentrazione= "<<dof_enum_C<<endl;
//ho fatto questo tipo di ciclo perché sul singolo elemento mi prendeva due volte il nodo iniziale: ora prendo entrambi nodi per il primo ciclo e poi per il secondo prendo solo il secondo valore del nodo
//NB: in questo modo metto in ordine la Cv per ogni branch: dato che la numerazione dei vertici salta per default di getfem, così ottengo la  numerazione corretta
		scalar_type sat=0;

		for(size_type j=0; j<dof_enum_H.size(); j++){
/*			cout<<"\n dof_enum_C["<<j<<"]= "<<dof_enum_C[j];
			cout<<"dof_enum_H["<<j<<"]= "<<dof_enum_H[j];
			cout<<"\n *************** "<<endl;
			cout<<"Cv[dof_enum_C["<<j<<"]]= "<<Cv[dof_enum_C[j]]<<endl;
			cout<<"Hi[dof_enum_H["<<j<<"]]= "<<Hi[dof_enum_H[j]]<<endl;
*/
			sat = dimensioning_saturation(Cv[dof_enum_C[j]]);
			psi[dof_enum_H[j]] = Hi[dof_enum_H[j]] * k1 * sat;
			//cout<<"psi[dof_enum_H["<<j<<"]]= "<<psi[dof_enum_H[j]]<<endl;
		}
/*
		vector_type saturation(mf_Hi[i].nb_dof()); gmm::clear(saturation);
		saturation = dimensioning_saturation(cv_i);

		for (size_type j=0; j<mf_Hi[i].nb_dof(); ++j)
		{
			//psi[j] = Hi[j]*k1*saturation[j];
			psi[j] = ht_i[j]*k1*saturation[j];
		}
*/
		//gmm::vscale(ht_i, saturation); //Ht*S(Cv);
		//gmm::scale(saturation, k1);	//Ht*k1*S(Cv) = PSI;

		getfem::interpolation(mf_Hi[i], mf_Uvi[i], psi, new_psi);

		pos=0;
	for(getfem::mr_visitor mrv(rg_branch); !mrv.finished(); ++mrv){
		//cout<<"\n Elemento della mesh fem della velocità per il calcolo di PSI: "<<mrv.cv()<<endl;
		if(pos==0){
			dof_enum_PSI.emplace_back(mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[0]); pos++;
			dof_enum_PSI.emplace_back(mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[1]); pos++;
			dof_enum_PSI.emplace_back(mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[2]); pos++;
			/*psi_i[pos] = new_psi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[0]];
			cout << "  psi_i[" << pos << "] = new_psi[" << mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[0] << "]" <<"="<<new_psi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[0]]<< endl;
			pos ++;
			psi_i[pos] = new_psi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[1]];
			cout << "  psi_i [" << pos << "] = new_psi[" << mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<new_psi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;
			psi_i[pos] = new_psi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[2]];
			cout << "  psi_i [" << pos << "] = new_psi[" << mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[2] << "]" <<"="<<new_psi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[2]]<< endl;
			pos ++;*/
			}
		else{
			dof_enum_PSI.emplace_back(mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[1]); pos++;
			dof_enum_PSI.emplace_back(mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[2]); pos++;
			/*psi_i[pos] = new_psi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[1]];
			cout << "  psi_i [" << pos << "] = new_psi[" << mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<new_psi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;
			psi_i[pos] = new_psi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[2]];
			cout << "  psi_i [" << pos << "] = new_psi[" << mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[2] << "]" <<"="<<new_psi[mf_Uvi[i].ind_basic_dof_of_element(mrv.cv())[2]]<< endl;
			pos ++;*/
			}
	}
//assemblo una nuova velocità: uv(1+PSI), ma PSI vive su mf_Hi (P1), invece la uv vive su mf_Uvi (P2) --> interpolo psi (che vive su mf_Hi[i]) sulla mesh
//della velocità mf_Uvi (P2), così non sto neanche a cambiare la mesh nella fase di assemblaggio

		for(size_type j=0; j<dof_enum_U.size(); j++){
			//cout<<"\n dof_enum_U["<<j<<"]= "<<dof_enum_U[j];
			//cout<<"dof_enum_PSI["<<j<<"]= "<<dof_enum_PSI[j];
			new_uvi[dof_enum_U[j]] = uvi[dof_enum_U[j]] * (1.0 + new_psi[dof_enum_PSI[j]]);
		}

		//cout<<"new_uvi= "<<new_uvi<<endl;
		return new_uvi;
	}; //end of modifing velocity


scalar_type 
oxygen_transport3d1d::dimensioning_saturation (scalar_type cv){
		/* la funzione assembla il vettore saturazione dimensionalizzato (dato che la cv del problema è adimensionale) per calcolare
			il vettore psi
		*/

		scalar_type coeff;
		coeff = param_oxy_transp.alpha_pl()/param_oxy_transp.C();
		scalar_type k2;
		k2 = pow(param_oxy_transp.Ps_50()*coeff,param_oxy_transp.delta()); //ml_O2/ml_B^delta

		scalar_type saturation;

		scalar_type num;
		num = pow(cv, param_oxy_transp.delta()-1.0);
		scalar_type denum;
		denum = pow(cv, param_oxy_transp.delta()) + k2;

		saturation = num/denum;

		return saturation;
	};// end of building saturation term
  
	  void oxygen_transport3d1d::assembly_oxy_transp (void)
	 {
 	 //1) Build the monolithic matrix AM
	 assembly_mat_oxy_transp();

	 //2) Build the monolithic rhs FM
	 assembly_rhs_oxy_transp();
	 }; // end of assembly
 
 
	void  
	oxygen_transport3d1d::assembly_mat_oxy_transp(void)
	{
	#ifdef M3D1D_VERBOSE_
	cout << "Allocating AM_oxy, UM_oxy, FM_oxy ..." << endl;
	#endif
	gmm::resize(AM_oxy, dof_oxy_transp.tot(), dof_oxy_transp.tot());	gmm::clear(AM_oxy);
	gmm::resize(UM_oxy, dof_oxy_transp.tot()); 		gmm::clear(UM_oxy);
	gmm::resize(FM_oxy, dof_oxy_transp.tot()); 		gmm::clear(FM_oxy);
	
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic matrix AM_oxy ..." << endl;
	#endif

	//Nel tessuto
	//Matrice di diffusione
	sparse_matrix_type Dt (dof_oxy_transp.Ct(), dof_oxy_transp.Ct()); gmm::clear(Dt);
	//Matrice di convezione
	sparse_matrix_type At (dof_oxy_transp.Ct(), dof_oxy_transp.Ct()); gmm::clear(At);
	//Matrice di reazione
	sparse_matrix_type Rt (dof_oxy_transp.Ct(), dof_oxy_transp.Ct()); gmm::clear(Rt);
	//Matrice per il drenaggio linfatico
	//sparse_matrix_type Lt (dof_oxy_transp.Ct(), dof_oxy_transp.Ct()); gmm::clear(Lt);

	//Nel vaso
	//Matrice di diffusione
	sparse_matrix_type Dv (dof_oxy_transp.Cv(), dof_oxy_transp.Cv()); gmm::clear(Dv);
	//Matrice di convezione
	sparse_matrix_type Av (dof_oxy_transp.Cv(), dof_oxy_transp.Cv()); gmm::clear(Av);

	//Matrici di scambio
	//Matrice di scambio tissue-to-tisse
	sparse_matrix_type Btt (dof_oxy_transp.Ct(), dof_oxy_transp.Ct()); gmm::clear(Btt);
	//Matrice di scambio vessel-to-tissue
	sparse_matrix_type Btv (dof_oxy_transp.Ct(), dof_oxy_transp.Cv()); gmm::clear(Btv);
	//Matrice di scambio tissue-to-vessel
	sparse_matrix_type Bvt (dof_oxy_transp.Cv(), dof_oxy_transp.Ct()); gmm::clear(Bvt);
	//Matrice di scambio vessel-to-vessel
	sparse_matrix_type Bvv (dof_oxy_transp.Cv(), dof_oxy_transp.Cv()); gmm::clear(Bvv);

	// Aux tissue-to-vessel averaging matrix
	sparse_matrix_type Mbar (dof_oxy_transp.Cv(), dof_oxy_transp.Ct()); gmm::clear(Mbar);

	// Aux tissue-to-vessel interpolation matrix
	sparse_matrix_type Mlin (dof_oxy_transp.Cv(), dof_oxy_transp.Ct()); gmm::clear(Mlin);
		
	#ifdef M3D1D_VERBOSE_
	cout<< "   Computing Peclet number..."<<endl;	
	#endif	
	//vectors containing the exact solution
	vector_type Ut(dof.Ut()); 	gmm::add(gmm::sub_vector(UM, gmm::sub_interval(0, dof.Ut())), 		Ut);
	vector_type Uv(dof.Uv()); 	gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())), Uv);
	//Interpolate Ut on polinomials of degree 0 
	pfem pf_U = fem_descriptor("FEM_PK(3,0)");
	mesh_fem mf_U(mesht);
	mf_U.set_qdim(3);
	mf_U.set_finite_element(mesht.convex_index(), pf_U);
	vector_type Ut_(mf_U.nb_dof());
	getfem::interpolation(mf_Ut, mf_U, Ut, Ut_);
	
	//RR
	//se c'è il trasporto di emoglobina (HEMOADVECTION) devo calcolare il nuovo Peclet con la velocità modificata
	//calcolo la velocità modificata, per ogni branch, per determinare peclet 
	if(descr_oxy_transp.HEMOADVECTION ==1){
		gmm::scale(Uv, 0.0);

		size_type shift = 0;
		size_type shift_h = 0;

		vector_type cv_guess(dof_oxy_transp.Cv(), param_oxy_transp.Cv_guess());
		for(size_type i; i<nb_branches; i++)
		{
			vector_type Hi(mf_Hi[i].nb_dof()); gmm::clear(Hi);
			vector_type Uvi(mf_Uvi[i].nb_dof()); gmm::clear(Uvi);
			vector_type new_uvi(mf_Uvi[i].nb_dof()); gmm::clear(new_uvi);

			if(i>0) shift += mf_Uvi[i-1].nb_dof();
			gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof())) ,  Uvi);

			if(i>0) shift_h += mf_Hi[i-1].nb_dof();
			gmm::add(gmm::sub_vector(UM_HT, gmm::sub_interval(shift_h, mf_Hi[i].nb_dof())), Hi);

			new_uvi = modifing_Uvi(Hi, Uvi, i, cv_guess);

			gmm::add(new_uvi , gmm::sub_vector(Uv, gmm::sub_interval(shift, mf_Uvi[i].nb_dof())));		
		}
	}

	scalar_type peclet_v= peclet(meshv, Uv, param_oxy_transp.Av(1), 1);
	scalar_type peclet_t= peclet(mesht, Ut_, param_oxy_transp.At(1), 3);
	
	#ifdef M3D1D_VERBOSE_
	cout<< "\nPeclet in vessels:    "<< peclet_v<<endl;
	cout<< "Peclet in tissue:    "<< peclet_t<<endl;
	#endif	

		
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Dt and Rt ..." << endl;
	#endif

/*
	//Coefficient for mass term:
	vector_type linf_coeff(dof.Pt()); gmm::clear(linf_coeff);
	vector_type Pl(dof.Pt(),PARAM.real_value("PL"));
	gmm::scale(Pl, -1.0); 
	gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut(), dof.Pt())) ,  linf_coeff);
	gmm::add(Pl ,  linf_coeff);
	gmm::scale (linf_coeff,param_transp.Q_pl(1));
*/

	//Coefficient for reaction term:
	vector_type consump_coeff(dof.Pt(), param_oxy_transp.M0());
	
	if(descr_oxy_transp.REACTION==1){
		cout<<" ******* Assembling Michaelis-Menten coefficient *******"<<endl;
	vector_type ct(dof.Pt(), param_oxy_transp.Ct_guess());
	vector_type PRESS50(dof.Pt(), param_oxy_transp.Pm_50());

	scalar_type k;
	k = param_oxy_transp.alpha_t()/param_oxy_transp.C();

	gmm::scale(PRESS50, k);

	gmm::add(PRESS50, ct);
	gmm::reciprocal(ct); //1/ct[i];
	gmm::scale(ct, param_oxy_transp.M0());
	gmm::copy(ct, consump_coeff);
	}
	
	//RR
	//Perché? cosi da avere la dimensione di cv_guess, consump_coeff e press50 identica a quella di UM_oxy(0, Ct) --> utile per il FPM, e 
	//dare come input a asm_tissue_transp il vettore di coefficienti (da costruire su mf_coeft) con le dimensioni corrette	

	//Build Dt, and Rt 
	asm_tissue_transp(Dt, Rt, mimt, mf_oxy_Ct, mf_coeft,  param_oxy_transp.At(), consump_coeff);

	// Check peclet number for instability
	if((descr_oxy_transp.ADVECTION==1) && (peclet_t>1))
		{ cout<<"WARNING!! Peclet > 1 in tissue: applying artificial diffusion"<<std::endl;	
	  	  gmm::scale(Dt, (1+peclet_t));}
	  
	
	gmm::add(Dt,
			 gmm::sub_matrix(AM_oxy, 
					gmm::sub_interval(0, dof_oxy_transp.Ct()), 
					gmm::sub_interval(0, dof_oxy_transp.Ct())));  
	
 	gmm::add(Rt, 
			  gmm::sub_matrix(AM_oxy, 
					gmm::sub_interval(0, dof_oxy_transp.Ct()), 
				 	gmm::sub_interval(0, dof_oxy_transp.Ct()))); 

/*		 	
	// Copy Lt: linfatic drainage in tissue
	 gmm::add(Lt, 
			  gmm::sub_matrix(AM_oxy, 
					gmm::sub_interval(0, dof_oxy_transp.Ct()), 
				 	gmm::sub_interval(0, dof_oxy_transp.Ct()))); 
*/

	
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Dv ..." << endl;
	#endif	
	
	// Build Dv
	asm_network_transp(Dv, mimv, mf_oxy_Cv, mf_coefv, param_oxy_transp.Av(), param.R());

	// Check peclet number for instability
	if((descr_oxy_transp.ADVECTION==1) && (peclet_v>1))
		{ cout<<"WARNING!! Peclet > 1 in network: applying artificial diffusion"<<endl;
   	 	  gmm::scale(Dv, (1+peclet_v));}
   	 	  cout<<"Artificial Diffusion = "<<(1+peclet_v)*param_oxy_transp.Av(1)<<endl;
	
	// Copy Dv: diffusion in network
	gmm::add(Dv, 	gmm::sub_matrix(AM_oxy, 
					gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv()), 
					gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())));
	
	if(descr_oxy_transp.ADVECTION ==0)
	{cout<<"No advection: only diffusion and reaction terms"<<endl;}
		
	else{		

	//ADVECTION	
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling At ..." << endl;
	#endif	
		
	//ADVECTION IN TISSUE		
	// Build At
	asm_advection_tissue(At, mimt, mf_oxy_Ct, mf_Ut, gmm::sub_vector(UM, gmm::sub_interval(0, dof.Ut())));

	// Copy At: advection in tissue
	gmm::add(At,
			  gmm::sub_matrix(AM_oxy, 
					gmm::sub_interval(0, dof_oxy_transp.Ct()), 
					gmm::sub_interval(0, dof_oxy_transp.Ct())));
	
 	#ifdef M3D1D_VERBOSE_
	cout<<"Assembling Av ..."<<endl;
	#endif


 	vector_type cv_guess (dof_oxy_transp.Cv(), param_oxy_transp.Cv_guess());

	size_type shift =0;
	size_type shift_h=0;
	
	for(size_type i=0; i<nb_branches; ++i)
	{
		vector_type Uvi(mf_Uvi[i].nb_dof()); gmm::clear(Uvi); //salvo la veclotià per ogni branch

		if(i>0) shift += mf_Uvi[i-1].nb_dof();
		gmm::copy(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof())) ,  Uvi);


		if(descr_oxy_transp.HEMOADVECTION==1){
			//cout<<"***********************"<<endl;
		vector_type Hi(mf_Hi[i].nb_dof()); gmm::clear(Hi); //salvo l'ematocrito per ogni branch (come è definito)

		if(i>0) shift_h += mf_Hi[i-1].nb_dof();
		gmm::copy(gmm::sub_vector(UM_HT, gmm::sub_interval(shift_h, mf_Hi[i].nb_dof())), Hi);

		vector_type new_uvi (mf_Uvi[i].nb_dof()); gmm::clear(new_uvi); //la nuova velocità: uv(1+PSI) --> quella da mettere nella fase di assemblaggio

		new_uvi = modifing_Uvi(Hi, Uvi, i, cv_guess);

		gmm::copy(new_uvi, Uvi);
		}

	asm_advection_network(Av, mimv, mf_oxy_Cv, mf_coefvi[i], mf_Uvi[i], mf_coefv,
							Uvi, param.lambdax(i), param.lambday(i), param.lambdaz(i),  param.R(), meshv.region(i) );
	}
	gmm::scale(Av, pi);

	// Copy Av: advection in network	
	gmm::add(Av,
			  gmm::sub_matrix(AM_oxy, 
					gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv()), 
					gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())));
	};// end of advection assembly


	bool COUPLING = PARAM.int_value("COUPLING", "flag for coupling-exchange term ");
	if(COUPLING==0)  { cout<< "Uncoupled problem: no exchange between tissue and vessels"<<endl; }
	else{
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling aux exchange matrices Mbar and Mlin ..." << endl;
	#endif

	if(PARAM.int_value("couple", "flag for coupling function (notaro 0, brambilla 1)")){
   		bool READ_INTERPOLATOR = PARAM.int_value("READ_INTERPOLATOR","flag for read interpolator from file"); //flag for read interpolator from file 
    	if (!READ_INTERPOLATOR){
		asm_exchange_aux_mat_transp(Mbar, Mlin, mimv, mf_oxy_Ct, mf_oxy_Cv, mf_coefv, param.R(), descr.NInt, nb_branches);
      	std::ostringstream mbar_name,mlin_name;
      	mbar_name << descr_oxy_transp.OUTPUT+"Mbar_transp";
      	mlin_name << descr_oxy_transp.OUTPUT+"Mlin_transp";
      	gmm::MatrixMarket_IO::write(mbar_name.str().c_str(), Mbar);
      	gmm::MatrixMarket_IO::write(mlin_name.str().c_str(), Mlin);
    	}
    	else
    	{
      	std::ostringstream mbar_name,mlin_name;
      	mbar_name << descr_oxy_transp.OUTPUT+"Mbar_transp";
      	mlin_name << descr_oxy_transp.OUTPUT+"Mlin_transp";
      	gmm::MatrixMarket_load(mbar_name.str().c_str(), Mbar);
      	gmm::MatrixMarket_load(mlin_name.str().c_str(), Mlin);
    	}
	}
	if(!PARAM.int_value("couple", "flag for coupling function (notaro 0, brambilla 1)")){
   		bool READ_INTERPOLATOR = PARAM.int_value("READ_INTERPOLATOR","flag for read interpolator from file");//flag for read interpolator from file 
    	if (!READ_INTERPOLATOR){

			asm_exchange_aux_mat(Mbar, Mlin, mimv, mf_oxy_Ct, mf_oxy_Cv, param.R(), descr.NInt);
      	std::ostringstream mbar_name,mlin_name;
      	mbar_name << "./vtk/Mbar";
      	mlin_name << "./vtk/Mlin";
      	gmm::MatrixMarket_IO::write(mbar_name.str().c_str(), Mbar);
      	gmm::MatrixMarket_IO::write(mlin_name.str().c_str(), Mlin);
    	}
    	else
    	{
      	std::ostringstream mbar_name,mlin_name;
      	mbar_name << "./vtk/Mbar";
      	mlin_name << "./vtk/Mlin";
      	gmm::MatrixMarket_load(mbar_name.str().c_str(), Mbar);
      	gmm::MatrixMarket_load(mlin_name.str().c_str(), Mlin);
    	}
	}
	
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling exchange matrices ..." << endl;
	#endif

	bool NEWFORM = PARAM.int_value("NEW_FORMULATION", "flag for the new formulation");
	
	// bluid oncotic term fot tissue
	vector_type ONCOTIC (dof.Pv());
	gmm::copy(gmm::sub_vector(UM, 
		gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())), ONCOTIC);
	gmm::mult_add(gmm::scaled(Mbar,-1.0), gmm::sub_vector(UM, 
		gmm::sub_interval(dof.Ut(), dof.Pt())),ONCOTIC);
		
	scalar_type picoef=param.sigma()*(param.pi_v()-param.pi_t());
        vector_type DeltaPi(dof.Pv(),picoef);

        gmm::add(gmm::scaled(DeltaPi,-1.0), ONCOTIC);	
	gmm::scale(ONCOTIC,0.5*(1.0-PARAM.real_value("sigma_oxy"))*param.Q(0)); //RR: aggiunta del coefficiente di riflessione soluto dipendente per l'O2

	if(descr_oxy_transp.TEST_ANALYTICAL){
		gmm::scale(ONCOTIC, 0.0) ;
	}
	
	// build permeability term for tissue
	vector_type PERM (dof.coefv());

	gmm::copy(param.R(), PERM);
	gmm::scale(PERM, 2*pi*param_oxy_transp.Y()[0]);

	//build exchange matrixes for tissue
	asm_exchange_mat_transp(Btt, Btv, Bvt, Bvv, mimv, mf_oxy_Cv, mf_coefv, mf_Pv, Mbar, Mlin, ONCOTIC, PERM, NEWFORM);

	// Copying Btt
	gmm::add(Btt,			 
			gmm::sub_matrix(AM_oxy, 
					gmm::sub_interval(0, dof_oxy_transp.Ct()), 
					gmm::sub_interval(0, dof_oxy_transp.Ct()))); 
	// Copying Btv

	gmm::add(Btv,
			gmm::sub_matrix(AM_oxy, 
					gmm::sub_interval(0, dof_oxy_transp.Ct()),
					gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv()))); 
	
	// Copying Bvt

	gmm::add(Bvt,								
			  gmm::sub_matrix(AM_oxy, 
			  		gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv()),
					gmm::sub_interval(0, dof_oxy_transp.Ct())));

	// Copying Bvv

	gmm::add(Bvv,
			  gmm::sub_matrix(AM_oxy, 
					gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv()), 
					gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())));	
	}	


	// De-allocate memory
	gmm::clear(Dt);    gmm::clear(Dv);
	gmm::clear(At);    gmm::clear(Av);
	gmm::clear(Rt);
	gmm::clear(Mbar);  gmm::clear(Mlin);
	gmm::clear(Btt);   gmm::clear(Btv);
	gmm::clear(Bvt);   gmm::clear(Bvv);

	}/* end of assembly_mat_transp */

	 

	void 
	oxygen_transport3d1d::assembly_rhs_oxy_transp(void)
	{
 
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic rhs FM_oxy ... " << endl;
	#endif
	
	#ifdef M3D1D_VERBOSE_
	cout << "  Building coupling dirichlet boundary term ..." << endl;
	#endif

	asm_coupled_bc_transp (AM_oxy, FM_oxy, mf_oxy_Ct, mf_oxy_Cv, BCt_oxy_transp, BCv_oxy_transp);
	
	//Impongo le BC per il tessuto	

	//Nel tessuto
	//Vettore al rhs per le BC nel tessuto
	vector_type Ft(dof_oxy_transp.Ct()); gmm::clear(Ft); //per le condizioni al contorno
	sparse_matrix_type Att(dof_oxy_transp.Ct(), dof_oxy_transp.Ct()); //per l'inetgrazione per parti del termine diffusivo

	gmm::add (gmm::sub_matrix(AM_oxy, 
			gmm::sub_interval(0,dof_oxy_transp.Ct()),
			gmm::sub_interval(0,dof_oxy_transp.Ct()))
			,Att);
	gmm::scale(	gmm::sub_matrix(AM_oxy, 
			gmm::sub_interval(0,dof_oxy_transp.Ct()),
			gmm::sub_interval(0,dof_oxy_transp.Ct()))
			,0.0);	
			
	gmm::add (gmm::sub_vector(FM_oxy, 
				gmm::sub_interval(0,dof_oxy_transp.Ct()))
			,Ft);

	gmm::scale(	 gmm::sub_vector(FM_oxy, 
				gmm::sub_interval(0,dof_oxy_transp.Ct()))
			,0.0);


	#ifdef M3D1D_VERBOSE_
	cout << "  Building tissue boundary term ..." << endl;
	#endif

	//Right Hand Side for tissue
	scalar_type beta_t  = PARAM.real_value("BETAtissue_transp", "Coefficient for mixed BC for transport problem in tissue");	
	asm_tissue_bc_transp(Ft, Att, mimt, mf_oxy_Ct, mf_coeft, BCt_oxy_transp, beta_t);
	//Att --> dall'integrazione per parti del termine diffusivo
	//Ft --> condizioni al contorno;

	gmm::add(Att, 
			gmm::sub_matrix(AM_oxy,
					gmm::sub_interval(0,dof_oxy_transp.Ct()),
					gmm::sub_interval(0,dof_oxy_transp.Ct())));
	gmm::add(Ft, 
			gmm::sub_vector(FM_oxy,
					gmm::sub_interval(0,dof_oxy_transp.Ct())));

	// De-allocate memory
	gmm::clear(Att);		
	gmm::clear(Ft);		


	//Impongo le BC per il vaso

	//Nel vaso
	//Vettore al rhs per le Bc nel vaso
	vector_type Fv (dof_oxy_transp.Cv()); gmm::clear(Fv);
	sparse_matrix_type Avv(dof_oxy_transp.Cv(), dof_oxy_transp.Cv());

	gmm::add(	gmm::sub_matrix(AM_oxy,
			gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv()),
			gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv()))
			, Avv);
	gmm::scale(	gmm::sub_matrix(AM_oxy,
			gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv()),
			gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv()))
			, 0.0);	
			
	gmm::add(gmm::sub_vector(FM_oxy,
				gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv()))
			,Fv);	 
	gmm::scale(	 gmm::sub_vector(FM_oxy,
				gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv()))
			,0.0);


	#ifdef M3D1D_VERBOSE_
	cout << "  Building vessel boundary term ..." << endl;
	#endif
		
	//Right Hand Side for vessels
	//Setting the BC for vessel
	scalar_type beta_v  = PARAM.real_value("BETAvessel_transp", "Coefficient for mixed BC for transport problem in vessels");
	asm_network_bc_transp(Fv, Avv, mimv, mf_oxy_Cv, mf_coefv, BCv_oxy_transp, beta_v, param.R());
	//Avv --> dall'integrazione per parti del termine diffusivo
	//Fv --> condizioni al contorno;

	gmm::add(Avv, 
			gmm::sub_matrix(AM_oxy,
					gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv()),
					gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv())));
	gmm::add(Fv, 
			gmm::sub_vector(FM_oxy,
					gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv())));
	// De-allocate memory
	gmm::clear(Avv);
	gmm::clear(Fv);	
	};// end of assembly_rhs_transp

	// Aux function for solver:
	// contains the list of different methods for solving (SuperLU, SAMG,GMRES, etc)
	bool oxygen_transport3d1d::solver (const size_type dof1, 
					   const size_type dof2,
					   const size_type dof3,
					   const size_type dof4){

	gmm::csc_matrix<scalar_type> A_transp;
	gmm::copy(AM_oxy, A_transp);
	
	vector_type F_transp(gmm::vect_size(FM_oxy));
	gmm::copy(FM_oxy, F_transp);


	if ( descr_oxy_transp.SOLVE_METHOD_OXY == "SuperLU" || descr_oxy_transp.SOLVE_METHOD_OXY == "SUPERLU" ) { // direct solver //
		#ifdef M3D1D_VERBOSE_
		cout << "  Applying the SuperLU method ... " << endl;
		#endif
		scalar_type cond;
		gmm::SuperLU_solve(A_transp, UM_oxy, F_transp, cond);
		cout << "  Condition number (transport problem): " << cond << endl;
	}
	else if(descr_oxy_transp.SOLVE_METHOD_OXY == "SAMG"){
	#ifdef WITH_SAMG	
	#ifdef M3D1D_VERBOSE_
		cout << "Solving the monolithic system ... " << endl;
	#endif

	//////////////////////////////////////AMG INTERFACE
	#ifdef M3D1D_VERBOSE_
	std::cout<<"converting A"<<std::endl;
	#endif
	gmm::csr_matrix<scalar_type> A_csr;


	int dim_matrix=dof1+dof2+dof3+dof4;
	gmm::copy(gmm::sub_matrix(AM_oxy,
			gmm::sub_interval(0 , dim_matrix),
			gmm::sub_interval(0 , dim_matrix)), A_csr);
	#ifdef M3D1D_VERBOSE_
	std::cout<<"converting X"<<std::endl;
	#endif
	std::vector<scalar_type> X,  B;

	gmm::resize(X,dim_matrix); gmm::clean(X, 1E-12);
	gmm::copy(gmm::sub_vector(UM_oxy,gmm::sub_interval(0,dim_matrix)),X);

	#ifdef M3D1D_VERBOSE_
	std::cout<<"converting B"<<std::endl;
	#endif
	gmm::resize(B,dim_matrix);gmm::clean(B, 1E-12);
	gmm::copy(gmm::sub_vector(FM_oxy,gmm::sub_interval(0,dim_matrix)),B);

	AMG amg("3d1d");
	amg.set_dof(dof1,dof2,dof3,dof4);

	amg.convert_matrix(A_csr);
	amg.solve(A_csr, X , B , 1);
	gmm::copy(amg.getsol(),gmm::sub_vector(UM_oxy,gmm::sub_interval(0,dim_matrix)));


	#ifdef SPARSE_INTERFACE
				for(int i = 0 ; i < nrows ; i++ ){U_1[i]=u[i];
					UM_transp[i]=u[i];	
				}
				gmm::copy(U_1, UM_oxy);
	#endif
	#ifdef CSC_INTERFACE
				for(int i = 0 ; i < nnu ; i++ ){
					U_2[i]=u_samg[i];UM_oxy[i]=u_samg[i];}
				gmm::copy(U_2,UM_oxy);
	#endif
	
	#else // with_samg=0
	std::cout<< "ERROR: you are trying to solve with samg, but WITH_SAMG=0"<<std::endl;
	std::cout<< "--> Call 'source configure.sh'and install the library again"<<std::endl;
	
	#endif
	}
	else { // Iterative solver //
/*
		// Iterations
		gmm::iteration iter(descr_oxy_transp.RES);  // iteration object with the max residu
		iter.set_noisy(1);               // output of iterations (2: sub-iteration)
		iter.set_maxiter(descr_oxy_transp.MAXITER); // maximum number of iterations

		// Preconditioners
		//! \todo Add preconditioner choice to param file
		// See \link http://download.gna.org/getfem/html/homepage/gmm/iter.html
		gmm::identity_matrix PM; // no precond
		//gmm::diagonal_precond<sparse_matrix_type> PM(AM); // diagonal preocond
		//gmm::ilu_precond<sparse_matrix_type> PM(AM);
		// ...
		//gmm::clear(AM);
		// See <http://download.gna.org/getfem/doc/gmmuser.pdf>, pag 15
	
		if ( descr_oxy_transp.SOLVE_METHOD == "CG" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Conjugate Gradient method ... " << endl;
			#endif
			gmm::identity_matrix PS;  // optional scalar product
			gmm::cg(AM_oxy, UM_oxy, F_oxy, PS, PM, iter);
		}
		else if ( descr_oxy_transp.SOLVE_METHOD == "BiCGstab" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the BiConjugate Gradient Stabilized method ... " << endl;
			#endif
			gmm::bicgstab(AM, UM, FM, PM, iter);
		}
		else if ( descr_oxy_transp.SOLVE_METHOD == "GMRES" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Generalized Minimum Residual method ... " << endl;
			#endif
			size_type restart = 50;
			gmm::gmres(AM_oxy, UM, FM, PM, restart, iter);
		}
		else if ( descr_oxy_transp.SOLVE_METHOD == "QMR" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Quasi-Minimal Residual method ... " << endl;
			#endif
			gmm::qmr(AM, UM, FM, PM, iter);
		}
		else if ( descr_oxy_transp.SOLVE_METHOD == "LSCG" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the unpreconditionned Least Square CG method ... " << endl;
			#endif
			gmm::least_squares_cg(AM, UM, FM, iter);
		}
		// Check
		if (iter.converged())
			cout << "  ... converged in " << iter.get_iteration() << " iterations." << endl;
		else if (iter.get_iteration() == descr_oxy_transp.MAXITER)
			cerr << "  ... reached the maximum number of iterations!" << endl;	*/

	}
		
	return true;
	} /* end of solver_oxy_transp */

////////////////////////////////////// SOLVER ////////////////////////////////////
	 bool oxygen_transport3d1d::solve_oxy_transp (void)
 	{
  	#ifdef M3D1D_VERBOSE_
	cout << "Solving the monolithic system ... " << endl;
	#endif
	
	
	double time = gmm::uclock_sec();
	//double time_partial = 0;
	//double time_count = 0;	
		

	// Solve the system on AM_oxy, UM_oxy, FM_oxy
	bool solved = solver(dof_oxy_transp.Ct(),dof_oxy_transp.Cv(),0,0);
	if (!solved) return false;
	

	//export solution
	#ifdef M3D1D_VERBOSE_
	std::cout<<"solved! going to export..."<<std::endl;
	#endif	
	
	if(descr_oxy_transp.REACTION && !descr_oxy_transp.HEMOADVECTION){
	string name_file = "_only_reaction";
	//std::ostringstream convert;
	//convert << time_count;
	//time_suff = convert.str();
	export_vtk_oxy_transp(name_file); 	
	}

	if(!descr_oxy_transp.REACTION && descr_oxy_transp.HEMOADVECTION){
	string name_file = "_only_psi";
	//std::ostringstream convert;
	//convert << time_count;
	//time_suff = convert.str();
	export_vtk_oxy_transp(name_file); 
	}

	if(descr_oxy_transp.REACTION && descr_oxy_transp.HEMOADVECTION){
	string name_file = "_3d1d";
	export_vtk_oxy_transp(name_file); 
	}
	if(descr_oxy_transp.TEST_ANALYTICAL){
		string name_file = "_analytical";
		export_vtk_oxy_transp(name_file);
	}

	/*
	string time_suff = "";
	std::ostringstream convert;
	convert << time_count;
	time_suff = convert.str();
	export_vtk_oxy_transp(); 
	*/

	#ifdef M3D1D_VERBOSE_		
	std::cout<<"exported!"<<std::endl;
	#endif	
	
	//if(!descr_oxy_transp.STATIONARY)
	//cout << "... time to solve : "	<< gmm::uclock_sec() - time_partial << " seconds\n";
	
/*	if(!descr_oxy_transp.STATIONARY){
	cout << endl<<"... time to solve all the time steps: " << gmm::uclock_sec() - time << " seconds\n";}
	else{
*/
	cout << endl<<"... time to solve : " << gmm::uclock_sec() - time << " seconds\n";
	
	return true;
 	}; // end of solve_oxy_transp


scalar_type 
computing_residual (vector_type Ct_new, vector_type Ct_old, vector_type Cv_new, vector_type Cv_old)
{

	gmm::scale(Ct_old, -1.0);
	gmm::add(Ct_old, Ct_new);

	gmm::scale(Cv_old, -1.0);
	gmm::add(Cv_old, Cv_new);

	scalar_type Tnew = gmm::vect_norm2(Ct_new);
	scalar_type Told = gmm::vect_norm2(Ct_old);

	scalar_type Vnew = gmm::vect_norm2(Cv_new);
	scalar_type Vold = gmm::vect_norm2(Cv_old);
	
	scalar_type error;
	error = Tnew/Told + Vnew/Vold;

	return error;
}


bool oxygen_transport3d1d::solve_oxygen_fixpoint (void)
{
	#ifdef  M3D1D_VERBOSE_
	cout<<"Start fix-point method..."<<endl;
	#endif

	bool RK;
	RK = solve_oxy_transp(); 

	clock_t time_G;

// Declaration of variables
	vector_type Ct_new(dof_oxy_transp.Ct()); gmm::clear(Ct_new);	//Ct(k)
	vector_type Ct_old(dof_oxy_transp.Ct()); gmm::clear(Ct_old);	//Ct(k-1)

	vector_type Cv_new(dof_oxy_transp.Cv()); gmm::clear(Cv_new);	//Cv(k)
	vector_type Cv_old(dof_oxy_transp.Cv()); gmm::clear(Cv_old);	//Cv(k-1)

	//Tissue matrices
	sparse_matrix_type At(dof_oxy_transp.Ct(), dof_oxy_transp.Ct()); gmm::clear(At);
	sparse_matrix_type Rt(dof_oxy_transp.Ct(), dof_oxy_transp.Ct()); gmm::clear(Rt);
	sparse_matrix_type Dt(dof_oxy_transp.Ct(), dof_oxy_transp.Ct()); gmm::clear(Dt);

	//Vessel Matrix
	sparse_matrix_type Av(dof_oxy_transp.Cv(), dof_oxy_transp.Cv()); gmm::clear(Av);
	sparse_matrix_type Dv(dof_oxy_transp.Cv(), dof_oxy_transp.Cv()); gmm::clear(Dv);

	//Exchange matrices
	sparse_matrix_type Btt (dof_oxy_transp.Ct(), dof_oxy_transp.Ct()); gmm::clear(Btt);
	sparse_matrix_type Btv (dof_oxy_transp.Ct(), dof_oxy_transp.Cv()); gmm::clear(Btv);
	sparse_matrix_type Bvt (dof_oxy_transp.Cv(), dof_oxy_transp.Ct()); gmm::clear(Bvt);
	sparse_matrix_type Bvv (dof_oxy_transp.Cv(), dof_oxy_transp.Cv()); gmm::clear(Bvv);

	sparse_matrix_type Mbar (dof_oxy_transp.Cv(), dof_oxy_transp.Ct()); gmm::clear(Mbar);
	sparse_matrix_type Mlin (dof_oxy_transp.Cv(), dof_oxy_transp.Ct()); gmm::clear(Mlin);


	//salvo la soluzione iniziale in Ct_old e in Cv_old
	gmm::copy(gmm::sub_vector(UM_oxy,
				gmm::sub_interval(0, dof_oxy_transp.Ct())), Ct_old);
	
	gmm::copy(gmm::sub_vector(UM_oxy,
				gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())), Cv_old);


	bool COUPLING = PARAM.int_value("COUPLING", "flag for coupling-exchange term ");
	if(COUPLING==0)  { cout<< "Uncoupled problem: no exchange between tissue and vessels"<<endl; }
	else{
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling aux exchange matrices Mbar and Mlin ..." << endl;
	#endif

	if(PARAM.int_value("couple", "flag for coupling function (notaro 0, brambilla 1)")){
   		bool READ_INTERPOLATOR = PARAM.int_value("READ_INTERPOLATOR","flag for read interpolator from file"); //flag for read interpolator from file 
    	if (!READ_INTERPOLATOR){
		asm_exchange_aux_mat_transp(Mbar, Mlin, mimv, mf_oxy_Ct, mf_oxy_Cv, mf_coefv, param.R(), descr.NInt, nb_branches);
      	std::ostringstream mbar_name,mlin_name;
      	mbar_name << descr_oxy_transp.OUTPUT+"Mbar_transp";
      	mlin_name << descr_oxy_transp.OUTPUT+"Mlin_transp";
      	gmm::MatrixMarket_IO::write(mbar_name.str().c_str(), Mbar);
      	gmm::MatrixMarket_IO::write(mlin_name.str().c_str(), Mlin);
    	}
    	else
    	{
      	std::ostringstream mbar_name,mlin_name;
      	mbar_name << descr_oxy_transp.OUTPUT+"Mbar_transp";
      	mlin_name << descr_oxy_transp.OUTPUT+"Mlin_transp";
      	gmm::MatrixMarket_load(mbar_name.str().c_str(), Mbar);
      	gmm::MatrixMarket_load(mlin_name.str().c_str(), Mlin);
    	}
	}
	if(!PARAM.int_value("couple", "flag for coupling function (notaro 0, brambilla 1)")){
   		bool READ_INTERPOLATOR = PARAM.int_value("READ_INTERPOLATOR","flag for read interpolator from file");//flag for read interpolator from file 
    	if (!READ_INTERPOLATOR){
			asm_exchange_aux_mat(Mbar, Mlin, mimv, mf_oxy_Ct, mf_oxy_Cv, param.R(), descr.NInt);
      	std::ostringstream mbar_name,mlin_name;
      	mbar_name << "./vtk/Mbar";
      	mlin_name << "./vtk/Mlin";
      	gmm::MatrixMarket_IO::write(mbar_name.str().c_str(), Mbar);
      	gmm::MatrixMarket_IO::write(mlin_name.str().c_str(), Mlin);
    	}
    	else
    	{
      	std::ostringstream mbar_name,mlin_name;
      	mbar_name << "./vtk/Mbar";
      	mlin_name << "./vtk/Mlin";
      	gmm::MatrixMarket_load(mbar_name.str().c_str(), Mbar);
      	gmm::MatrixMarket_load(mlin_name.str().c_str(), Mlin);
    	}
	}
}// end of coupling

	//compute peclet
	vector_type Ut(dof.Ut()); 	gmm::add(gmm::sub_vector(UM, gmm::sub_interval(0, dof.Ut())), 		Ut);
	
	//Interpolate Ut on polinomials of degree 0 
	pfem pf_U = fem_descriptor("FEM_PK(3,0)");
	mesh_fem mf_U(mesht);
	mf_U.set_qdim(3);
	mf_U.set_finite_element(mesht.convex_index(), pf_U);
	vector_type Ut_(mf_U.nb_dof());
	getfem::interpolation(mf_Ut, mf_U, Ut, Ut_);

	//computing peclet
	scalar_type peclet_t= peclet(mesht, Ut_, param_oxy_transp.At(1), 3);

	//Impongo un residuo (per la concentrazion) e un massimo di iterazioni 
	scalar_type oxyres = descr_oxy_transp.Residual_OXY;
	scalar_type err = 1.0; 

	//Impongo un residuo (per il bilanzio di massa)
	//scalar_type toll_diff = ...; es. 1.0e-3;
	//scalar_type toll_adv = ...; es. 1.0e-7
	//come condizione per il while --> count_diff==Jv_oxy_transp.size() && count_adv==Jv_oxy_transp.size()

	int max_iter = descr_oxy_transp.Max_iterations_OXY;
	int iteration=1;

	//da fare solo nel caso di non linearità
	if(!descr_oxy_transp.HEMOADVECTION && !descr_oxy_transp.REACTION){
		cout<<"Any non-linearities are in the system"<<endl;
	}

	else {
		time_G=clock();
while (RK && iteration<=max_iter && err>oxyres)
{

	if(descr_oxy_transp.REACTION==1){
		
 	#ifdef M3D1D_VERBOSE_
	cout<<"Assembling Rt in FPM..."<<endl;
	#endif

	gmm::scale(gmm::sub_matrix(AM_oxy, 
							gmm::sub_interval(0, dof_oxy_transp.Ct()),
							gmm::sub_interval(0, dof_oxy_transp.Ct())),0.0);
	vector_type ct_guess(dof_oxy_transp.Ct());
	vector_type consump_coeff(dof.Pt());

	//Riga di aggiornamento
	gmm::copy(Ct_old, ct_guess);

	getfem::interpolation(mf_oxy_Ct, mf_coeft, ct_guess, consump_coeff);

	scalar_type k;
	k = param_oxy_transp.alpha_t()/param_oxy_transp.C();

	vector_type PRESS50(dof.Pt(), param_oxy_transp.Pm_50());

	gmm::scale(PRESS50, k);

	gmm::add(PRESS50, consump_coeff);
	gmm::reciprocal(consump_coeff);
	gmm::scale(consump_coeff, param_oxy_transp.M0());

	asm_tissue_transp(Dt, Rt, mimt, mf_oxy_Ct, mf_coeft,  param_oxy_transp.At(), consump_coeff);

	if((descr_oxy_transp.ADVECTION==1) && (peclet_t>1))
		{ cout<<"WARNING!! Peclet > 1 in tissue: applying artificial diffusion"<<std::endl;	
	  	  gmm::scale(Dt, (1+peclet_t));}

	asm_advection_tissue(At, mimt, mf_oxy_Ct, mf_Ut, gmm::sub_vector(UM, gmm::sub_interval(0, dof.Ut())));

	if(COUPLING==1){
	bool NEWFORM = PARAM.int_value("NEW_FORMULATION", "flag for the new formulation");
	// bluid oncotic term fot tissue
	vector_type ONCOTIC (dof.Pv());
	gmm::copy(gmm::sub_vector(UM, 
		  		  gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())),
		  ONCOTIC);
	gmm::mult_add(gmm::scaled(Mbar,-1.0), 
		  gmm::sub_vector(UM, 
		  		  gmm::sub_interval(dof.Ut(), dof.Pt())),
		  ONCOTIC);
	scalar_type picoef=param.sigma()*(param.pi_v()-param.pi_t());
        vector_type DeltaPi(dof.Pv(),picoef);
        gmm::add(gmm::scaled(DeltaPi,-1.0), ONCOTIC);	
	gmm::scale(ONCOTIC,0.5*(1.0-PARAM.real_value("sigma_oxy"))*param.Q(0));
	
	// build permeability term for tissue
	vector_type PERM (dof.coefv());
	gmm::copy(param.R(), PERM);
	gmm::scale(PERM, 2*pi*param_oxy_transp.Y()[0]);

	asm_exchange_mat_transp(Btt, Btv, Bvt, Bvv, mimv, mf_oxy_Cv, mf_coefv, mf_Pv, Mbar, Mlin, ONCOTIC, PERM, NEWFORM);

	gmm::add(Btt, gmm::sub_matrix(AM_oxy, 
									gmm::sub_interval(0, dof_oxy_transp.Ct()), 
									gmm::sub_interval(0, dof_oxy_transp.Ct())));
	}

		gmm::add(Dt, gmm::sub_matrix(AM_oxy,
									gmm::sub_interval(0,dof_oxy_transp.Ct()),
									gmm::sub_interval(0,dof_oxy_transp.Ct())));

		gmm::add(Rt, gmm::sub_matrix(AM_oxy,
									gmm::sub_interval(0,dof_oxy_transp.Ct()),
									gmm::sub_interval(0,dof_oxy_transp.Ct())));

		gmm::add(At, gmm::sub_matrix(AM_oxy,
									gmm::sub_interval(0,dof_oxy_transp.Ct()),
									gmm::sub_interval(0,dof_oxy_transp.Ct())));

		gmm::clear(At);
		gmm::clear(Dt);
		gmm::clear(Rt);
		gmm::clear(Btt); gmm::clear(Bvv);
		gmm::clear(Btv); gmm::clear(Bvt);

		gmm::clear(gmm::sub_vector(FM_oxy, gmm::sub_interval(0, dof_oxy_transp.Ct())));
 	}

if(descr_oxy_transp.HEMOADVECTION==1){

 	#ifdef M3D1D_VERBOSE_
	cout<<"Assembling Av in FPM..."<<endl;
	#endif

	vector_type Uv(dof.Uv()); gmm::clear(Uv);

	gmm::scale( gmm::sub_matrix(AM_oxy,
					gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv()),
					gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())), 0.0);

 	vector_type cv_guess (dof_oxy_transp.Cv()); gmm::clear(cv_guess);
 	gmm::copy(Cv_old, cv_guess); //Riga di aggiornamento

	size_type shift =0;
	size_type shift_h=0;

	for(size_type i=0; i<nb_branches; ++i)
	{
		vector_type Hi(mf_Hi[i].nb_dof()); gmm::clear(Hi);
		vector_type Uvi(mf_Uvi[i].nb_dof()); gmm::clear(Uvi);
		vector_type new_uvi (mf_Uvi[i].nb_dof()); gmm::clear(new_uvi); 

		if(i>0) shift_h += mf_Hi[i-1].nb_dof();
		if(i>0) shift += mf_Uvi[i-1].nb_dof();

		gmm::copy(gmm::sub_vector(UM_HT, 
			gmm::sub_interval(shift_h, mf_Hi[i].nb_dof())), Hi);

		gmm::copy(gmm::sub_vector(UM,
			gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof())), Uvi);

		new_uvi = modifing_Uvi(Hi, Uvi, i, cv_guess); //new velocity per ogni branch i-esimo

	gmm::add(new_uvi, gmm::sub_vector(Uv, gmm::sub_interval(shift, mf_Uvi[i].nb_dof()))); //salvo la velocità modificata

	asm_advection_network(Av, mimv, mf_oxy_Cv, mf_coefvi[i], mf_Uvi[i], mf_coefv, 
							new_uvi, param.lambdax(i), param.lambday(i), param.lambdaz(i),  param.R(), meshv.region(i));
	}
	gmm::scale(Av, pi);
	
	scalar_type peclet_v = peclet(meshv, Uv, param_oxy_transp.Av(1), 1);

	asm_network_transp(Dv, mimv, mf_oxy_Cv, mf_coefv, param_oxy_transp.Av(), param.R());

	if((descr_oxy_transp.ADVECTION==1) && (peclet_v>1))
		{ cout<<"WARNING!! Peclet > 1 in network: applying artificial diffusion= "<<1+peclet_v<<std::endl;	
	  	  gmm::scale(Dv, (1+peclet_v));}

	if(COUPLING==1){
	bool NEWFORM = PARAM.int_value("NEW_FORMULATION", "flag for the new formulation");
	// bluid oncotic term fot tissue
	vector_type ONCOTIC (dof.Pv());
	gmm::copy(gmm::sub_vector(UM, 
		  		  gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())),
		  ONCOTIC);
	gmm::mult_add(gmm::scaled(Mbar,-1.0), 
		  gmm::sub_vector(UM, 
		  		  gmm::sub_interval(dof.Ut(), dof.Pt())),
		  ONCOTIC);
	scalar_type picoef=param.sigma()*(param.pi_v()-param.pi_t());
        vector_type DeltaPi(dof.Pv(),picoef);
        gmm::add(gmm::scaled(DeltaPi,-1.0), ONCOTIC);	
	gmm::scale(ONCOTIC,0.5*(1.0-PARAM.real_value("sigma_oxy"))*param.Q(0));

	vector_type PERM (dof.coefv());
	gmm::copy(param.R(), PERM);
	gmm::scale(PERM, 2*pi*param_oxy_transp.Y()[0]);

	asm_exchange_mat_transp(Btt, Btv, Bvt, Bvv, mimv, mf_oxy_Cv, mf_coefv, mf_Pv, Mbar, Mlin, ONCOTIC, PERM, NEWFORM);

	gmm::add(Bvv, gmm::sub_matrix(AM_oxy, 
			gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv()),
			gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())));	
	}

	gmm::add(Dv, gmm::sub_matrix(AM_oxy, 
			gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv()),
			gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())));		

	gmm::add(Av, gmm::sub_matrix(AM_oxy, 
			gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv()),
			gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())));	

	gmm::clear(Av);
	gmm::clear(Dv);
	gmm::clear(Bvv); gmm::clear(Btt);
	gmm::clear(Btv); gmm::clear(Bvt);

	gmm::clear(gmm::sub_vector(FM_oxy, gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())));
	} //fine HEMOADVECTION

		//Assemblo nuovamente il rhs
		//gmm::clear(FM_oxy);
		assembly_rhs_oxy_transp();
		//cout<<"Risolvo il sistema nuovamente"<<endl;
			RK = solve_oxy_transp();

			gmm::copy(gmm::sub_vector(UM_oxy,
							gmm::sub_interval(0, dof_oxy_transp.Ct())), Ct_new);

			gmm::copy(gmm::sub_vector(UM_oxy,
							gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())), Cv_new);
		//calcolare l'errore tra il passo k e quello k-1
			err = computing_residual (Ct_new, Ct_old, Cv_new, Cv_old);

	#ifdef M3D1D_VERBOSE_
	cout<<"Updating the solution"<<endl;
	#endif
	//under-relaxation
/*	
	if(param_oxy_transp.underOXY()!=1){
	gmm::scale(Cv_new,param_oxy_transp.underOXY());
	gmm::scale(Cv_old,(1-param_oxy_transp.underOXY()));
	gmm::add(Cv_old,Cv_new);
	gmm::scale(Cv_old,1/(1-param_oxy_transp.underOXY())); //?
	}
*/
	gmm::copy(Ct_new, Ct_old);
	gmm::copy(Cv_new, Cv_old);

	#ifdef M3D1D_VERBOSE_
	cout << "Checking Residuals at iteration: "<< iteration << "..." << endl;
	#endif
	//Saving residual values in an output file
	//SaveResidual << iteration << "\t" << resSol << "\t" << resCM << "\t" << resH << endl;

			if(descr_oxy_transp.PRINT_RESIDUALS)  {
			cout << "\nStep n°:" << iteration << "\nSolution Residual = " << err <<endl;
			//cout << "\t\t\t\tTime: " <<  ((float)t)/CLOCKS_PER_SEC << " s "<< endl;
					}
			cout << "********************************************************" << endl;

	iteration++;
	gmm::clear(Ct_new);
	gmm::clear(Cv_new);
	}//esco dal while


	time_G=clock()-time_G;
	cout<< "Iterative Process Time = " << ((float)time_G)/CLOCKS_PER_SEC << " s"<< endl;
} //fine dell'if (presenza o no di non linearità)

	return true;
}; //end of oxygen fixpoint
	

	//Compute the residuals for mass balance at each junction 
	void oxygen_transport3d1d::mass_balance(void){

	#ifdef M3D1D_VERBOSE_		
	cout << " Compute MBD and MBA "   << endl;
	#endif	

	//PROVA
	/*
	cout<<"Dimensione Jv_oxy_transp= "<<Jv_oxy_transp.size()<<endl;
	// initialize the MBD and MBA to zero (clear at eac time step)
	for (size_type i=0; i<Jv_oxy_transp.size(); ++i){
		Jv_oxy_transp[i].MBD=0;
		Jv_oxy_transp[i].MBA=0;

		//PROVA
		cout<<"Analizziamo cosa c'è dentro a Jv_oxy_transp"<<endl;
		cout<<"\n Jv_oxy_transp["<<i<<"].label= "<<Jv_oxy_transp[i].label<<endl;
		cout<<"\n Jv_oxy_transp["<<i<<"].rg= "<<Jv_oxy_transp[i].rg<<endl;
		cout<<"\n Jv_oxy_transp["<<i<<"].idx= "<<Jv_oxy_transp[i].idx<<endl;
		cout<<"\n Jv_oxy_transp["<<i<<"].branches= "<<Jv_oxy_transp[i].branches<<endl;
		///////
	}	
	*/
	size_type shift = 0; //counter for branches
	
	for (size_type i=0; i<mf_Uvi.size(); ++i){ // branch loop 	
		if(i>0) shift += mf_Uvi[i-1].nb_dof(); 
		mesh_region &rg_branch = meshv.region(i); // branch region

		for (size_type j=0; j<Jv_oxy_transp.size(); ++j){ //junction loop
			mesh_region &rg_junction = meshv.region(Jv_oxy_transp[j].rg);  // junction region: dico qual è la regione della giunzione j-esima
			// Iterators for all the branches which flow in the junction j
			std::vector<long signed int>::const_iterator bb = Jv_oxy_transp[j].branches.begin();
			std::vector<long signed int>::const_iterator be = Jv_oxy_transp[j].branches.end();
		
			//Check if outflow of branch i is in junction j: per quale branch la regione della giunzione j-esima è l'outflow		
			if ((std::find(bb, be, +i) != be)){

				// find last element of the branch
				getfem::mr_visitor ii(rg_branch); int temp=0;
				for(; !ii.finished() && temp==0; ++ii){ // loop for convexes of the branch i
					getfem::mr_visitor ii_temp=ii;
					++ii_temp;
					if(ii_temp.finished()) temp=1;
				} 
				
				//import radius of the branch
				scalar_type Ri = param.R(mimv, i);
				// find the dof for last point for Cv e Uv
				size_type last_C, last_C2, last_U;
				vector_type dof_enum_C,dof_enum_U;
				int fine_C=0, fine_U=0;
				for (mr_visitor mrv(rg_branch); !mrv.finished(); ++mrv){
				for (auto b : mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv()))
					{dof_enum_C.emplace_back(b);
					fine_C++;}			
				for (auto b : mf_Uvi[i].ind_basic_dof_of_element(mrv.cv()))
					{dof_enum_U.emplace_back(b);
					fine_U++;}			
				}

				last_C = dof_enum_C[fine_C-1];
				last_C2= dof_enum_C[fine_C-2];
				last_U = dof_enum_U[fine_U-1];
				dof_enum_C.clear(); dof_enum_U.clear();
	
				scalar_type DIFF=0, ADV=0;

				//Compute the diffusive flux
				DIFF = pi* Ri * Ri*(
						UM_oxy[dof_oxy_transp.Ct()+last_C]-UM_oxy[dof_oxy_transp.Ct()+last_C2] )
						/estimate_h(meshv, ii.cv()) ;
				//Compute the advective fluxes
				ADV = pi* Ri * Ri*UM[dof.Ut()+dof.Pt()+shift+last_U]*UM_oxy[dof_oxy_transp.Ct()+last_C];
	
				#ifdef M3D1D_VERBOSE_		
				cout << "------------------------------------------ "   << endl;
				cout <<"in branch "<< i << " and junction "<< j <<" (region number :  "<< Jv_oxy_transp[j].rg<<" )"<<endl;
				cout<<"MBD_partial = "<< DIFF<< endl;
				cout<<"MBA_partial = "<< ADV<<  endl;
				#endif	
		
				Jv_oxy_transp[j].MBD -= DIFF;
				Jv_oxy_transp[j].MBA -= ADV;
					
			}// end of check for outflow branch
	
			//Check if inflow of branch i is in junction j			
			if ( (i!=0) &&  (std::find(bb, be, -i) != be )){  //notice that, in build_vessel_transp, we assume that the branch zero cannot have the inflow in a junction. 


				// find first element of the branch
				getfem::mr_visitor ii(rg_branch);
			
				//import radius of the branch
				scalar_type Ri = param.R(mimv, i);
				// find the dof for last point for Cv e Uv
				size_type first_C, first_C2, first_U;
				vector_type dof_enum_C,dof_enum_U;
				int fine_C=0, fine_U=0;
				for (mr_visitor mrv(rg_branch); !mrv.finished(); ++mrv){
				for (auto b : mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv()))
					{dof_enum_C.emplace_back(b);
					fine_C++;}			
				for (auto b : mf_Uvi[i].ind_basic_dof_of_element(mrv.cv()))
					{dof_enum_U.emplace_back(b);
					fine_U++;}			
				}
	
				first_C = dof_enum_C[0];
				first_C2= dof_enum_C[1];
				first_U = dof_enum_U[0];
				dof_enum_C.clear(); dof_enum_U.clear();
	
				scalar_type DIFF=0,ADV=0;
	
				//Compute the diffusive flux: Fick's Law
				 DIFF = pi* Ri * Ri*(
						UM_oxy[dof_oxy_transp.Ct()+first_C2]-UM_oxy[dof_oxy_transp.Ct()+first_C] )
						/estimate_h(meshv, ii.cv()) ;
	
	
				//Compute the advective fluxes: mass flow rate
				ADV = pi* Ri * Ri*UM[dof.Ut()+dof.Pt()+shift+first_U]*UM_oxy[dof_oxy_transp.Ct()+first_C];
	
				#ifdef M3D1D_VERBOSE_		
				cout << "------------------------------------------ "   << endl;
				cout <<"in branch "<< i << " and junction "<< j <<" (region number :  "<< Jv_oxy_transp[j].rg<<" )"<<endl;
				cout<<"MBD_partial = "<< DIFF<< endl;
				cout<<"MBA_partial = "<< ADV<<  endl;
				#endif	
		
				Jv_oxy_transp[j].MBD += DIFF;
				Jv_oxy_transp[j].MBA += ADV;
					
			}// end of check for outflow branch

		} //end of junction loop
	} // end of branch loop	


/* DA RIMETTERE
	cout << "  Junctions: " << endl;
	for (size_type i=0; i<Jv_oxy_transp.size(); ++i){
		cout << "    -  label=" << Jv_oxy_transp[i].label 
			 << ", value=" << Jv_oxy_transp[i].value << ", ind=" << Jv_oxy_transp[i].idx 
			 << ", rg=" << Jv_oxy_transp[i].rg << ", branches=" << Jv_oxy_transp[i].branches << endl; 
		cout << " Mass balance of diffusive fluxes = " << Jv_oxy_transp[i].MBD << endl; 
		cout << " Mass balance of advective fluxes = " << Jv_oxy_transp[i].MBA << endl;
		cout << "             ------------------- "   << endl;
	} 	
		cout << "----------------------------------------------- "   << endl;
*/

	
	}; // end of mass_balance

vector_type
oxygen_transport3d1d::computing_oxyhemoglobin(vector_type Hi, size_type i, vector_type Cv){
//cout<<"\n ********** Calcolando PSI ..."<<endl;
//cout<<"Branch #"<<i<<endl;
	scalar_type k1;
		k1 = param_oxy_transp.MCHC() * param_oxy_transp.N() / param_oxy_transp.C();
	scalar_type k2;
		k2 = pow(param_oxy_transp.Ps_50()*param_oxy_transp.alpha_pl(),param_oxy_transp.delta())/pow(param_oxy_transp.C(),param_oxy_transp.delta()); //ml_O2/ml_B^delta

		vector_type dof_enum_C, dof_enum_H;

		vector_type psi(Hi.size()); gmm::clear(psi);
		

		mesh_region &rg_branch = meshv.region(i);	

		size_type pos=0;
	for (getfem::mr_visitor mrv(rg_branch); !mrv.finished(); ++mrv)
		{
		if(pos == 0)
			{
			dof_enum_C.emplace_back(mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0]); pos++;
			dof_enum_C.emplace_back(mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]); pos++;
			/*cv_i[pos] = Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0]];
			cout << " cv_i [" << pos << "] = cv[" << mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0] << "]" <<"="<<Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0]]<< endl;
			pos ++;
			cv_i[pos] = Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " cv_i [" << pos << "] = cv[" << mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;*/
			}
		else{
			dof_enum_C.emplace_back(mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]); pos++;
			/*cv_i[pos] = Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " cv_i [" << pos << "] = cv[" << mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;*/
			}
		}

	pos=0;
	for(getfem::mr_visitor mrv(rg_branch); !mrv.finished(); ++mrv){
		//cout<<"\n Elemento della mesh fem dell'ematocrito: "<<mrv.cv()<<endl;
		if(pos==0){
			dof_enum_H.emplace_back(mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0]); pos++;
			dof_enum_H.emplace_back(mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]); pos++;
			/*ht_i[pos] = Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0]];
			cout << " ht_i [" << pos << "] = Hi[" << mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0] << "]" <<"="<<Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0]]<< endl;
			pos ++;
			ht_i[pos] = Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " ht_i [" << pos << "] = Hi[" << mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;*/
			}
			else{
			dof_enum_H.emplace_back(mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]); pos++;
			/*ht_i[pos] = Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " ht_i [" << pos << "] = Hi[" << mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;*/
			}
	}

	scalar_type sat;
	scalar_type num;
	scalar_type denum;
	for(size_type j=0; j<dof_enum_H.size(); j++){
		num = pow(Cv[dof_enum_C[j]], param_oxy_transp.delta());
		denum = num + k2;
		sat = num/denum;
		psi[dof_enum_H[j]] = Hi[dof_enum_H[j]] * k1 * sat;
	}

/*
		//gmm::reciprocal(denum);
		//gmm::vscale(num, denum); //calcolo della saturazione: OUTPUT: denum

		for (size_type j=0; j<Hi.size(); j++)
		{
			psi[j] = ht_i[j]*k1*saturation[j];
			cout<<"Psi["<<j<<"]= "<<psi[j]<<endl;
		}

		pos=0;
	for(getfem::mr_visitor mrv(rg_branch); !mrv.finished(); ++mrv){
		cout<<"\n Elemento della mesh per il calcolo di PSI: "<<mrv.cv()<<endl;
		if(pos==0){
			psi_i[pos] = psi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0]];
			cout << " psi_i [" << pos << "] = Psi[" << mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0] << "]" <<"="<<psi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0]]<< endl;
			pos ++;
			psi_i[pos] = psi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " psi_i [" << pos << "] = Psi[" << mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<psi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;
			}
			else{
			psi_i[pos] = psi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " psi_i [" << pos << "] = Psi[" << mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<psi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;
			}
	}
	*/
		//gmm::vscale(ht_i, denum);
		//gmm::scale(denum, k1);
		//gmm::copy(denum, psi);

	return psi;
}

///////////////////////////////////////////////////////////////////////////////////////////////
////////////////		CALCOLO DI CVTOT 		///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
vector_type
oxygen_transport3d1d::computing_totaloxygen(vector_type Hi, size_type i, vector_type Cv){
//cout<<"\n ********** Calcolando CVTOT ..."<<endl;
//cout<<"Branch #"<<i<<endl;
	scalar_type k1;
		k1 = param_oxy_transp.MCHC() * param_oxy_transp.N() / param_oxy_transp.C();
	scalar_type k2;
		k2 = pow(param_oxy_transp.Ps_50()*param_oxy_transp.alpha_pl(),param_oxy_transp.delta())/pow(param_oxy_transp.C(),param_oxy_transp.delta()); //ml_O2/ml_B^delta

		vector_type dof_enum_C, dof_enum_H;

		vector_type psi(Hi.size()); gmm::clear(psi);

		vector_type Cvtot(Hi.size()); gmm::clear(Cvtot);
		
		/*
		vector_type cv_i(Hi.size()); gmm::clear(cv_i);	
		vector_type ht_i(Hi.size()); gmm::clear(ht_i);
		vector_type psi_i(Hi.size()); gmm::clear(psi_i);
		*/

		mesh_region &rg_branch = meshv.region(i);	

		size_type pos=0;
	for (getfem::mr_visitor mrv(rg_branch); !mrv.finished(); ++mrv)
		{
		if(pos == 0)
			{
			dof_enum_C.emplace_back(mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0]); pos++;
			dof_enum_C.emplace_back(mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]); pos++;
			/*cv_i[pos] = Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0]];
			cout << " cv_i [" << pos << "] = cv[" << mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0] << "]" <<"="<<Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0]]<< endl;
			pos ++;
			cv_i[pos] = Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " cv_i [" << pos << "] = cv[" << mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;*/
			}
		else{
			dof_enum_C.emplace_back(mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]); pos++;
			/*cv_i[pos] = Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " cv_i [" << pos << "] = cv[" << mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<Cv[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;*/
			}
		}

	pos=0;
	for(getfem::mr_visitor mrv(rg_branch); !mrv.finished(); ++mrv){
		//cout<<"\n Elemento della mesh fem dell'ematocrito: "<<mrv.cv()<<endl;
		if(pos==0){
			dof_enum_H.emplace_back(mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0]); pos++;
			dof_enum_H.emplace_back(mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]); pos++;
			/*ht_i[pos] = Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0]];
			cout << " ht_i [" << pos << "] = Hi[" << mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0] << "]" <<"="<<Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0]]<< endl;
			pos ++;
			ht_i[pos] = Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " ht_i [" << pos << "] = Hi[" << mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;*/
			}
			else{
			dof_enum_H.emplace_back(mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]); pos++;
			/*ht_i[pos] = Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " ht_i [" << pos << "] = Hi[" << mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<Hi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;*/
			}
	}

	scalar_type sat;
	scalar_type num;
	scalar_type denum;
	for(size_type j=0; j<dof_enum_H.size(); j++){
		num = pow(Cv[dof_enum_C[j]], param_oxy_transp.delta());
		denum = num + k2;
		sat = num/denum;
		psi[dof_enum_H[j]] = Hi[dof_enum_H[j]] * k1 * sat;

		Cvtot[dof_enum_H[j]] = Cv[dof_enum_C[j]] + psi[dof_enum_H[j]];
	}

/*
		//gmm::reciprocal(denum);
		//gmm::vscale(num, denum); //calcolo della saturazione: OUTPUT: denum

		for (size_type j=0; j<Hi.size(); j++)
		{
			psi[j] = ht_i[j]*k1*saturation[j];
			cout<<"Psi["<<j<<"]= "<<psi[j]<<endl;
		}

		pos=0;
	for(getfem::mr_visitor mrv(rg_branch); !mrv.finished(); ++mrv){
		cout<<"\n Elemento della mesh per il calcolo di PSI: "<<mrv.cv()<<endl;
		if(pos==0){
			psi_i[pos] = psi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0]];
			cout << " psi_i [" << pos << "] = Psi[" << mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0] << "]" <<"="<<psi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[0]]<< endl;
			pos ++;
			psi_i[pos] = psi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " psi_i [" << pos << "] = Psi[" << mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<psi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;
			}
			else{
			psi_i[pos] = psi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " psi_i [" << pos << "] = Psi[" << mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<psi[mf_Hi[i].ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;
			}
	}
	*/
		//gmm::vscale(ht_i, denum);
		//gmm::scale(denum, k1);
		//gmm::copy(denum, psi);

	//return psi;
		return Cvtot;
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
 

 const void oxygen_transport3d1d::export_vtk_oxy_transp (const string & time_suff,const string & suff)
 {

  if (PARAM.int_value("VTK_EXPORT"))
  {
	#ifdef M3D1D_VERBOSE_
	cout << "Exporting the solution (vtk format) to " << descr.OUTPUT << " ..." << endl;
	#endif
	#ifdef M3D1D_VERBOSE_
	cout << "  Saving the results from the monolithic unknown vector ... " << endl;
	#endif
	
	// Array of unknown dof of the interstitial velocity
	vector_type Ct(dof_oxy_transp.Ct()); 

	// Array of unknown dof of the network velocity
	vector_type Cv(dof_oxy_transp.Cv()); 

	vector_type Cv_dim(dof_oxy_transp.Cv());
	
	//Copy solution
	gmm::copy(gmm::sub_vector(UM_oxy, 
		gmm::sub_interval(0, dof_oxy_transp.Ct())),  Ct);
	gmm::copy(gmm::sub_vector(UM_oxy, 
		gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())), Cv);	

	gmm::copy(Cv, Cv_dim);

 	gmm::scale(Cv, param_oxy_transp.C());

 	gmm::scale(Ct, param_oxy_transp.C());

/* DA RIMETTERE
	for(size_type a=0;a<dof_oxy_transp.Ct();a++){
		cout<<"Ct["<<a<<"]= "<<Ct[a]<<endl;
	}

	for(size_type b=0;b<dof_oxy_transp.Cv();b++){
		cout<<"Cv["<<b<<"]= "<<Cv[b]<<endl;
	}
*/

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Ct ..." << endl;
	#endif
	vtk_export exp_Ct(descr_oxy_transp.OUTPUT+"Ct"+time_suff+".vtk");
	exp_Ct.exporting(mf_oxy_Ct);
	exp_Ct.write_mesh();
	exp_Ct.write_point_data(mf_oxy_Ct, Ct, "Ct");


	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Cv ..." << endl;
	#endif
	vtk_export exp_Cv(descr_oxy_transp.OUTPUT+"Cv"+time_suff+".vtk");
	exp_Cv.exporting(mf_oxy_Cv);
	exp_Cv.write_mesh();
	exp_Cv.write_point_data(mf_oxy_Cv, Cv, "Cv");

	cout<<"NUMERI DOF DELLA MESH T E V"<<endl;
	cout<<"#DOF di T: "<<dof_oxy_transp.Ct()<<endl;
	cout<<"#DOF di V: "<<dof_oxy_transp.Cv()<<endl;

/* 	DA RIMETTERE
	if(descr_oxy_transp.HEMOADVECTION==1){
	#ifdef M3D1D_VERBOSE_
	cout<< "  Exporting Oxyhemoglobin and Total Oxygen Concentration.." <<endl;
	#endif
	size_type start = 0;
	size_type length = 0;
	size_type k, last_dof;
	last_dof=nb_branches*mf_Hi[nb_branches-1].nb_dof();
	for(size_type i=0; i<nb_branches; ++i){
		vector_type psi(mf_Hi[i].nb_dof());
		vector_type Cvtot(mf_Hi[i].nb_dof());
		vector_type Hi(mf_Hi[i].nb_dof());

		if(i>0) start += mf_Hi[i-1].nb_dof();
		length = mf_Hi[i].nb_dof();

		gmm::copy(gmm::sub_vector(UM_HT, gmm::sub_interval(start, length)), Hi);
		psi = computing_oxyhemoglobin(Hi, i, Cv_dim);
		Cvtot = computing_totaloxygen(Hi, i, Cv_dim);
		gmm::scale(psi, param_oxy_transp.C());
		gmm::scale(Cvtot, param_oxy_transp.C());
		//cout<<"psi= "<<psi<<endl;
		//cout<<"Cvtot= "<<Cvtot<<endl;

		vtk_export exp_Oh(descr_oxy_transp.OUTPUT+"Oh"+suff+std::to_string(i)+".vtk");
			exp_Oh.exporting(mf_Hi[i]);
			exp_Oh.write_mesh();
			exp_Oh.write_point_data(mf_Hi[i], psi, "Psi"); 

		vtk_export exp_Cvtot(descr_oxy_transp.OUTPUT+"Cvtot"+suff+std::to_string(i)+".vtk");
			exp_Cvtot.exporting(mf_Hi[i]);
			exp_Cvtot.write_mesh();
			exp_Cvtot.write_point_data(mf_Hi[i], Cvtot, "Cvtot"); 
			}
	}
	*/

	#ifdef M3D1D_VERBOSE_
	cout << "... export done, visualize the data file with (for example) Paraview" << endl; 
	#endif
  }
 }; // end of export_transp

/*


  // Interface with problem3d1d class
  	//! Initialize the problem
	void transport3d1d::init_fluid (int argc, char *argv[])
	{ 
	problem3d1d::init(argc, argv);
	};
	
	//! Assemble the problem
	void transport3d1d::assembly_fluid (void)
	{ 
	problem3d1d::assembly();
	};
	//! Solve the problem
	bool transport3d1d::solve_fluid (void)
	{ 
	return problem3d1d::solve();
	};	
	//! Export the solution
	const void transport3d1d::export_vtk_fluid (const string & suff)
	{ 
	problem3d1d::export_vtk(suff);
	};
*/
 } // end of namespace
