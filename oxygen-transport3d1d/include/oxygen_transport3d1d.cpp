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
  @date   September 2016 - September 2018.
  @brief Definizione della classe principale per risolvere il problema del trasporto di ossigeno
 */
 
 #include <oxygen_transport3d1d.hpp>
 #include <AMG_Interface.hpp>
 #include <cmath>
 #include "gmm/gmm_inoutput.h"
 #include "getfem/getfem_import.h"


 namespace getfem {


	// Inizializzo il problema
 	void oxygen_transport3d1d::init_oxy_transp (int argc, char *argv[]) 
 	{
 	#ifdef M3D1D_VERBOSE_
 	std::cout << "initialize transport problem..."<<std::endl<<std::endl;
 	#endif

	//PARAM.read_command_line(argc, argv);
	//1. Import data (algorithm specifications, boundary conditions, ...)	
	//import_data_oxy_transp();	//OK!
	//2. Import mesh for tissue (3D) and vessel network (1D)
	build_mesh_oxy_transp();	//OK!	
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


	mf_Hi.clear();
	mf_coefh.clear();

	problemHT::set_im_and_fem();

	}; // end of set_im_and_oxy_fem
	
	
	// Build problem parameters
	void oxygen_transport3d1d::build_param_oxy_transp(void)
	{
	
	#ifdef M3D1D_VERBOSE_
	cout << "Building parameters for tissue and vessel problems ..." << endl;
	#endif
	param.build(PARAM, mf_coeft, mf_coefv,mf_coefvi); 
	param_oxy_transp.build(PARAM, mf_coeft, mf_coefv);
	//param_oxy_transp.build_oxy(PARAM, mf_coeft_oxy);

	//PROVA::
	cout<<"Dof di mf_coeft: "<<mf_coeft.nb_dof()<<endl;
	cout<<"Dof di mf_coefv: "<<mf_coefv.nb_dof()<<endl;
	cout<<"Dof di Pv: "<<dof.Pv()<<endl;
	cout<<"Dof di Pt: "<<dof.Pt()<<endl;
	cout<<"Dof di Cv"<<dof_oxy_transp.Cv()<<endl;
	cout<<"Dof di Ct"<<dof_oxy_transp.Ct()<<endl;
	///////////////


	#ifdef M3D1D_VERBOSE_
	cout << param_oxy_transp ;
	#endif

	cout<<"C'è il test analtico? TEST_ANALYTICAL = "<<descr_oxy_transp.TEST_ANALYTICAL<<endl;
	cout<<"Coefficiente di diffusione test = "<<param_oxy_transp.Dt()<<endl;


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
for(size_type j=0; j<5; j++){
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
      cout << "Convex of index " << i << endl;
      bgeot::pconvex_structure cvs = mesht.structure_of_convex(i);
      cout << "Number of vertices: " << cvs->nb_points() << endl;
      cout << "Number of faces: " << cvs->nb_faces() << endl;
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

/*
		PROVA: lista dei convessi di meshv e rispettivi indici dei nodi
		cout<<"convesso numero: "<<cv<<endl;
		cout<<"i0= "<<i0<<endl;
		cout<<"i1= "<<i1<<endl;
*/
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

	} /* end of build_vessel_boundary_oxy_transp */

	//////////////////////////// ASSEMBLAGGIO ////////////////////////////
  
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
	//Matrcei di convezione
	sparse_matrix_type Av (dof_oxy_transp.Cv(), dof_oxy_transp.Cv()); gmm::clear(Av);

	//Matrici di scambio
	//Matrice di scambio tissue-to-tisse
	sparse_matrix_type Btt (dof_oxy_transp.Ct(), dof_oxy_transp.Ct()); gmm::clear(Btt);
	//Matrice di scambio tissue-to-vessel
	sparse_matrix_type Btv (dof_oxy_transp.Ct(), dof_oxy_transp.Cv()); gmm::clear(Btv);
	//Matrice di scambio vessel-to-tissue
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
	//compute peclet
	scalar_type peclet_v= peclet(meshv, Uv, param_oxy_transp.Av(1), 1);
	scalar_type peclet_t= peclet(mesht, Ut_, param_oxy_transp.At(1), 3);
	
	#ifdef M3D1D_VERBOSE_
	cout<< "Peclet in vessels:    "<< peclet_v<<endl;
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
	vector_type ct_guess(dof_oxy_transp.Ct(), param_oxy_transp.Ct_guess());
	
	vector_type PRESS50(dof_oxy_transp.Ct(), PARAM.real_value("Pm_50"));
	gmm::scale(PRESS50, param_oxy_transp.alpha_t());
	
	gmm::add(PRESS50, ct_guess);

	for (size_type i=0; i<dof_oxy_transp.Ct(); i++)
	{
	ct_guess[i] = param_oxy_transp.M0()/ct_guess[i];
	}

	//interpolo il coefficiente di consumo (ct_guess) su un polinomio di grado 0
	pfem pf_consump = fem_descriptor("FEM_PK(3,0)");
	mesh_fem mf_R(mesht);
	mf_R.set_qdim(1);
	mf_R.set_finite_element(mesht.convex_index(), pf_consump);
	vector_type consump_coeff(mf_R.nb_dof());

	getfem::interpolation(mf_oxy_Ct, mf_R, ct_guess, consump_coeff);

	cout<<"Dimensione mf_R: "<<mf_R.nb_dof()<<endl;
	cout<<"Dimensione mf_coeft: "<<mf_coeft.nb_dof()<<endl;
	cout<<"consump_coeff = "<<consump_coeff[0]<<endl;

	//RR
	//getfem::interpolation(mf_oxy_Ct, mf_coeft, ct_guess, consump_coeff);
	//CONSUMP_COEFF = 0.3099;
	//Perché? cosi da avere la dimensione di cv_guess, consump_coeff e press50 identica a quella di UM_oxy(0, Ct) --> utile per il FPM, e 
	//dare come input a asm_tissue_transp il vettore di coefficienti (da costruire su mf_coeft) con le dimensioni corrette	

	//Build Dt, and Rt 
	asm_tissue_transp(Dt, Rt, mimt, mf_oxy_Ct, mf_coeft,  param_oxy_transp.At(), consump_coeff);
	  
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
	cout<<"***************************Raggio: "<<param.R()<<endl;

	// Check peclet number for instability
	 if((descr_oxy_transp.ADVECTION==1) && (peclet_v>1))
		{ cout<<"WARNING!! Peclet > 1 in network: applying artificial diffusion"<<endl;
   	 	  gmm::scale(Dv, (1+peclet_v)); }
		
	// Copy Dv: diffusion in network 	
	gmm::add(Dv, 	gmm::sub_matrix(AM_oxy, 
					gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv()), 
					gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())));
	


	if(descr_oxy_transp.ADVECTION ==0)
	{cout<<"No advection: only diffusion and reaction terms"<<endl;}
		
	else{		

	//ADVECTION	
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling At and Av ..." << endl;
	#endif	
	
		
	//ADVECTION IN TISSUE		
	// Build At
	asm_advection_tissue(At, mimt, mf_oxy_Ct, mf_Ut, descr_oxy_transp.TEST_ANALYTICAL, gmm::sub_vector(UM, gmm::sub_interval(0, dof.Ut())));

	// Copy At: advection in tissue
	gmm::add(At,
			  gmm::sub_matrix(AM_oxy, 
					gmm::sub_interval(0, dof_oxy_transp.Ct()), 
					gmm::sub_interval(0, dof_oxy_transp.Ct()))); 	
		
	// build Av
	size_type shift = 0;
	for(size_type i=0; i<nb_branches; ++i){

		if(i>0) shift += mf_Uvi[i-1].nb_dof();

		vector_type Uvi (mf_Uvi[i].nb_dof()); gmm::clear(Uvi);
		gmm::add(gmm::sub_vector(UM, 
				gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof())) ,  Uvi);


		asm_advection_network(Av, mimv, mf_oxy_Cv, mf_coefvi[i], mf_Uvi[i], mf_coefv, Uvi, param.lambdax(i), param.lambday(i), param.lambdaz(i),  param.R(), meshv.region(i) );
	}
	gmm::scale(Av, pi);

	// Copy Av: advection in network
	gmm::add(Av,
			  gmm::sub_matrix(AM_oxy, 
					gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv()), 
					gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())));
	}/* end of advection assembly*/ 


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
		  		  gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())),
		  ONCOTIC);
		
	gmm::mult_add(gmm::scaled(Mbar,-1.0), 
		  gmm::sub_vector(UM, 
		  		  gmm::sub_interval(dof.Ut(), dof.Pt())),
		  ONCOTIC);
		
	scalar_type picoef=param.sigma()*(param.pi_v()-param.pi_t());
        vector_type DeltaPi(dof.Pv(),picoef);
        gmm::add(gmm::scaled(DeltaPi,-1.0), ONCOTIC);	
	gmm::scale(ONCOTIC,0.5*(1.0-param.sigma())*param.Q(0));

	if(descr_oxy_transp.TEST_ANALYTICAL){
		gmm::scale(ONCOTIC, 0.0);
	}
	
	// build permeability term for tissue
	vector_type PERM (dof.coefv());
	gmm::copy(param.R(), PERM);
	gmm::scale(PERM, 2*pi*param_oxy_transp.Y()[0]);
	// PERM adimensionale: 27.5 (Perm/U = 27.5e-3/1.0e-4)

	//build exchange matrixes for tissue
	asm_exchange_mat_transp(Btt, Btv, Bvt, Bvv,
			mimv, mf_oxy_Cv, mf_coefv, mf_Pv, Mbar, Mlin, 
			ONCOTIC, PERM, NEWFORM);


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


	cout<<"FM_oxy size: "<<FM_oxy.size()<<endl;
	cout<<"AM_oxy rows: "<<gmm::mat_nrows(AM_oxy)<<" e AM_oxy cols: "<<gmm::mat_ncols(AM_oxy)<<endl;
	cout<<"Dof di mf_oxy_Ct: "<<mf_oxy_Ct.nb_dof()<<endl;
	cout<<"Dof di mf_oxy_Cv: "<<mf_oxy_Cv.nb_dof()<<endl;

	//PROVA:
	for(size_type bc; bc<BCt_oxy_transp.size();bc++){
	cout<<"BCt_transp LABEL= "<<BCt_oxy_transp[bc].label<<endl;
	cout<<"BCt_transp VALUE= "<<BCt_oxy_transp[bc].value<<endl;
	//cout<<"BCt_transp IDX= "<<BCt_oxy_transp[bc].idx<<endl;
	cout<<"BCt_transp REGION= "<<BCt_oxy_transp[bc].rg<<"\n"<<endl;
	}
	/*
	cout<<"******************************************************"<<endl;
	for(size_type bc; bc<BCv_oxy_transp.size();bc++){
	cout<<"BCv_transp LABEL= "<<BCv_oxy_transp[bc].label<<endl;
	cout<<"BCv_transp VALUE= "<<BCv_oxy_transp[bc].value<<endl;
	//cout<<"BCv_transp IDX= "<<BCv_oxy_transp[bc].idx<<endl;
	cout<<"BCv_transp REGION= "<<BCv_oxy_transp[bc].rg<<"\n"<<endl;
	}*/
	/////////////////

	asm_coupled_bc_transp (AM_oxy, FM_oxy, mf_oxy_Ct, mf_oxy_Cv, BCt_oxy_transp, BCv_oxy_transp);
	
	
	//Impongo le BC per il tessuto	

	//Nel tessuto
	//Vettore al rhs per le BC nel tessuto
	vector_type Ft(dof_oxy_transp.Ct()); gmm::clear(Ft); //per le condizioni al contorno
	sparse_matrix_type Att(dof_oxy_transp.Ct(), dof_oxy_transp.Ct()); //per l'inetgrazione per parti del termine diffusivo


//	vector_type FM_temp(dof_oxy_transp.tot()); gmm::clear(FM_temp);
//	sparse_matrix_type AM_temp(dof_oxy_transp.tot(), dof_oxy_transp.tot()); gmm::clear(AM_temp);


//	gmm::copy(AM_oxy, AM_temp);
//	gmm::copy(FM_oxy, FM_temp);


	gmm::add (gmm::sub_matrix(AM_oxy, //AM_temp
			gmm::sub_interval(0,dof_oxy_transp.Ct()),
			gmm::sub_interval(0,dof_oxy_transp.Ct()))
			,Att);
	gmm::scale(	gmm::sub_matrix(AM_oxy, //AM_temp
			gmm::sub_interval(0,dof_oxy_transp.Ct()),
			gmm::sub_interval(0,dof_oxy_transp.Ct()))
			,0.0);	
			
	gmm::add (gmm::sub_vector(FM_oxy, //FM_temp
				gmm::sub_interval(0,dof_oxy_transp.Ct()))
			,Ft);	
	gmm::scale(	 gmm::sub_vector(FM_oxy, //FM_temp
				gmm::sub_interval(0,dof_oxy_transp.Ct()))
			,0.0);

	#ifdef M3D1D_VERBOSE_
	cout << "  Building tissue boundary term ..." << endl;
	#endif
	
	//Right Hand Side for tissue
	scalar_type beta_t  = PARAM.real_value("BETAtissue_transp", "Coefficient for mixed BC for transport problem in tissue");	
	
	asm_tissue_bc_transp(Ft, Att, mimt, mf_oxy_Ct, mf_coeft, BCt_oxy_transp, beta_t);
	//cout<<"Ho impostato le BC (DIR e MIX) per il tessuto"<<endl;
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

	gmm::add(	gmm::sub_matrix(AM_oxy, //AM_temp
			gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv()),
			gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv()))
			, Avv);
	gmm::scale(	gmm::sub_matrix(AM_oxy, //AM_temp
			gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv()),
			gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv()))
			, 0.0);	
			
	gmm::add(	 gmm::sub_vector(FM_oxy, //AM_temp
				gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv()))
			,Fv);	 
	gmm::scale(	 gmm::sub_vector(FM_oxy, //FM_temp
				gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv()))
			,0.0);


	#ifdef M3D1D_VERBOSE_
	cout << "  Building vessel boundary term ..." << endl;
	#endif
		
	//Right Hand Side for vessels
	//Setting the BC for vessel
	scalar_type beta_v  = PARAM.real_value("BETAvessel_transp", "Coefficient for mixed BC for transport problem in vessels");
	asm_network_bc_transp(Fv, Avv, mimv, mf_oxy_Cv, mf_coefv, BCv_oxy_transp, beta_v, param.R());
	cout<<"Ho impostato le BC (DIR e MIX) per il vaso"<<endl;
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
	


	/*

	//Assembling Ov: oxyhemoglobin advcetion vector
	//Vettore per il trasporto l'ossiemoglobina (termine non lineare)
	vector_type Ov (dof_oxy_transp.Cv()); gmm::clear(Ov);	

	vector_type cv_guess (dof_oxy_transp.Cv(), param_oxy_transp.Cv_guess());
					
	size_type shift =0;
	size_type shift_h=0;
	
	scalar_type k1;
	k1= param_oxy_transp.N()*param_oxy_transp.MCHC();	
	scalar_type k2;
	k2 = pow((param_oxy_transp.Ps_50()*param_oxy_transp.alpha_t()),param_oxy_transp.delta());		
	for(size_type i=0; i<nb_branches; ++i){

		cout<<"Sono nel ciclo!!!!!"<<endl;
	cout<<"mf_Hi size: "<<mf_Hi[i].nb_dof() <<endl;
	cout<<"mf_oxy_Cv size: "<<mf_oxy_Cv.nb_dof()<<endl;

		vector_type cv_i(mf_Hi[i].nb_dof()); gmm::clear(cv_i);
		vector_type Hi(mf_Hi[i].nb_dof()); gmm::clear(Hi);
		vector_type Uvi(mf_Uvi[i].nb_dof()); gmm::clear(Uvi);
		vector_type psi(mf_Hi[i].nb_dof()); gmm::clear(psi);
		vector_type power(mf_Hi[i].nb_dof()); gmm::clear(power);


		/*size_type pos=0;
		for (getfem::mr_visitor mrv(mf_oxy_Cv.linked_mesh().region(i)); !mrv.finished(); ++mrv)
		for (auto b : mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv()))
			{
	cout<<"Sono con ind_basic_dof_of_element!!!!!"<<endl;
	cout<<"b= "<<b<<endl;
			cv_i[pos] = cv_guess[b]; //ottengo cv_i --> la concentrazione definita ramo per ramo
			pos++;
			}

	size_type pos=0;
	for (getfem::mr_visitor mrv(mf_oxy_Cv.linked_mesh().region(i)); !mrv.finished(); ++mrv){
		if(pos == 0){
			cv_i[pos] = cv_guess[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0]];
			cout << " cv_i [" << pos << "] = cv_guess [" << mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0] << "]" << endl;
			pos ++;
			cv_i[pos] = cv_guess[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " cv_i [" << pos << "] = cv_guess [" << mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1] << "]" << endl;
			pos ++;
			}
		else{
			cv_i[pos] = cv_guess[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]];
			cout << " cv_i [" << pos << "] = cv_guess [" << mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1] << "]" << endl;
			pos ++;

			}
	//ho fatto questo tipo di ciclo perché sul singolo elemento mi prendeva due volte il nodo iniziale: ora prendo entrambi nodi per il primo ciclo e poi per il secondo prendo solo il secondo valore del nodo
	}

		if(i>0) shift_h += mf_Hi[i-1].nb_dof();
		if(i>0) shift += mf_Uvi[i-1].nb_dof();
		
		gmm::add(gmm::sub_vector(UM_HT, 
			gmm::sub_interval(shift_h, mf_Hi[i].nb_dof())), Hi);
		gmm::add(gmm::sub_vector(UM,
			gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof())) ,  Uvi);

		cout<<"psi size: "<<psi.size()<<endl;
		cout<<"cv_i size: "<<cv_i.size()<<endl;

		for(size_type f; f<cv_i.size(); f++){
		cout<<"valore di cv_i: "<<cv_i[f]<<endl;
		}

		for (size_type j=0; j<mf_Hi[i].nb_dof(); j++)
		{
			cout<<"cv_i["<<j<<"]="<<cv_i[j]<<endl;
			psi[j] = Hi[j]*k1*pow(cv_i[j], param_oxy_transp.delta())/(pow(cv_i[j], param_oxy_transp.delta())+k2);
			
		}
	asm_hemoadvection_rhs_network(Ov, mimv, mf_oxy_Cv, mf_coefvi[i], mf_Uvi[i], mf_coefv, mf_Hi[i], Uvi, param.lambdax(i), param.lambday(i), param.lambdaz(i),  param.R(), psi, meshv.region(i));	
	}
	
	gmm::copy(Ov, gmm::sub_vector(FM_oxy, 
					gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv())));
	*/


	// De-allocate memory
	//gmm::clear(AM_temp);
	//gmm::clear(FM_temp);
	}/* end of assembly_rhs_transp */

	 
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

	cout<<"Residuo di Ct= "<<Tnew/Told<<endl;
	cout<<"Residuo di Cv= "<<Vnew/Vold<<endl;
	
	scalar_type error;

	error = Tnew/Told + Vnew/Vold;

	return error;
}


bool oxygen_transport3d1d::solve_oxygen_fixpoint (void)
{
	#ifdef  M3D1D_VERBOSE_
	cout<<"Start FixPoint method..."<<endl;
	#endif

	cout<<"Risolvo il sistema"<<endl;
	bool RK;
	RK = solve_oxy_transp(); 

	vector_type Ct(dof_oxy_transp.Ct()); 
	vector_type Cv(dof_oxy_transp.Cv()); 


	//Copy solution
	gmm::copy(gmm::sub_vector(UM_oxy, 
		gmm::sub_interval(0, dof_oxy_transp.Ct())),  Ct);
	gmm::copy(gmm::sub_vector(UM_oxy, 
		gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())), Cv);
	cout<<"VALORI DI Cv e Ct SENZA FPM"<<endl;
	for(size_type a; a<dof_oxy_transp.Ct();a++)
		cout<<"Ct["<<a<<"]= "<<Ct[a]<<endl;

	for(size_type b; b<dof_oxy_transp.Cv();b++)
		cout<<"Cv["<<b<<"]= "<<Cv[b]<<endl;



// 2 - declaration of variables

	//vettori per la concentrazione

	vector_type Ct_new(dof_oxy_transp.Ct()); gmm::clear(Ct_new);	//Ct(k)
	vector_type Ct_old(dof_oxy_transp.Ct()); gmm::clear(Ct_old);	//Ct(k-1)

	vector_type Cv_new(dof_oxy_transp.Cv()); gmm::clear(Cv_new);	//Cv(k)
	vector_type Cv_old(dof_oxy_transp.Cv()); gmm::clear(Cv_old);	//Cv(k-1)

	//Definisco le matrici che non cambiano durante il while
	sparse_matrix_type Rt_temp (dof_oxy_transp.Ct(), dof_oxy_transp.Ct());
	//sparse_matrix_type Dt (dof_oxy_transp.Ct(), dof_oxy_transp.Ct());
	sparse_matrix_type Dt_temp (dof_oxy_transp.Ct(), dof_oxy_transp.Ct());

	//sparse_matrix_type At (dof_oxy_transp.Ct(), dof_oxy_transp.Ct());
	sparse_matrix_type At_temp (dof_oxy_transp.Ct(), dof_oxy_transp.Ct());

	//sparse_matrix_type Btt (dof_oxy_transp.Ct(), dof_oxy_transp.Ct());
	sparse_matrix_type Btt_temp (dof_oxy_transp.Ct(), dof_oxy_transp.Ct());

	sparse_matrix_type Btv (dof_oxy_transp.Ct(), dof_oxy_transp.Cv());
	sparse_matrix_type Bvt (dof_oxy_transp.Cv(), dof_oxy_transp.Ct());
	sparse_matrix_type Bvv (dof_oxy_transp.Cv(), dof_oxy_transp.Cv());
	sparse_matrix_type Mbar (dof_oxy_transp.Cv(), dof_oxy_transp.Ct());
	sparse_matrix_type Mlin (dof_oxy_transp.Cv(), dof_oxy_transp.Ct());

/*
	vector_type Ft (dof_oxy_transp.Ct());
	vector_type Fv (dof_oxy_transp.Cv());

	gmm::copy(
			gmm::sub_vector(FM_oxy,
				gmm::sub_interval(0, dof_oxy_transp.Ct())), Ft);

	gmm::copy(
			gmm::sub_vector(FM_oxy,
				gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())), Fv);
*/

	//Inizializzo la marice Rt e il vettore Ov che cambino nel while
	sparse_matrix_type Rt (dof_oxy_transp.Ct(), dof_oxy_transp.Ct());
	vector_type Ov(dof_oxy_transp.Cv());

	//definisco il polinomio di grado 0 su cui interpolare il coefficiente di consumo di ossigeno (in Rt)
	pfem pf_consump = fem_descriptor("FEM_PK(3,0)");
	mesh_fem mf_R(mesht);
	mf_R.set_qdim(1);
	mf_R.set_finite_element(mesht.convex_index(), pf_consump);
	
	//Assemblo Dt
	getfem::asm_stiffness_matrix_for_laplacian(Dt_temp ,mimt,mf_oxy_Ct, mf_coeft, param_oxy_transp.At());//, mesh_region::all_convexes());
	
	
	//lo aggiungo all'interno del while

	//Assemblo At
	asm_advection_tissue(At_temp, mimt, mf_oxy_Ct, mf_Ut, descr_oxy_transp.TEST_ANALYTICAL, gmm::sub_vector(UM, gmm::sub_interval(0, dof.Ut())));

	//lo aggiungo all'interno del while

	//Assemblo le matrici di scambio Btt, Btv, Bvt e Bvv (sarebbe meglio solo Btt)
	bool COUPLING = PARAM.int_value("COUPLING", "flag for coupling-exchange term ");
	if(COUPLING==0)  { cout<< "Uncoupled problem: no exchange between tissue and vessels"<<endl; }
	else{

	if(PARAM.int_value("couple", "flag for coupling function (notaro 0, brambilla 1)")){
   		bool READ_INTERPOLATOR = PARAM.int_value("READ_INTERPOLATOR", "flag for read interpolator from file ");
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
   		bool READ_INTERPOLATOR = PARAM.int_value("READ_INTERPOLATOR", "flag for read interpolator from file ");
    	if (!READ_INTERPOLATOR){
    		//con couple=0 e read_interpolator=0;
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
	gmm::scale(ONCOTIC,(1.0-param.sigma())*param.Q(0));

	if(descr_oxy_transp.TEST_ANALYTICAL){
		gmm::scale(ONCOTIC,0.0);
	}
	
	// build permeability term for tissue
	vector_type PERM (dof.coefv());
	gmm::copy(param.R(), PERM);
	gmm::scale(PERM, 2*pi*param_oxy_transp.Y()[0]);


	//build exchange matrixes for tissue
	asm_exchange_mat_transp(Btt_temp, Btv, Bvt, Bvv,
			mimv, mf_oxy_Cv, mf_coefv, mf_Pv, Mbar, Mlin, 
			ONCOTIC, PERM, NEWFORM);

	
	//lo aggiunego all'interno del while
	}

	//Copio i blocchi di AM_oxy e FM_oxy nelle matrici e nei vettori che non cambino nel while
	/*
	gmm::copy(gmm::sub_matrix(AM_oxy,
						gmm::sub_interval(0, dof_oxy_transp.Ct()),
						gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())),Btv);

	gmm::copy(gmm::sub_matrix(AM_oxy,
						gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv()),
						gmm::sub_interval(0, dof_oxy_transp.Ct())),Bvt);

	gmm::copy(gmm::sub_matrix(AM_oxy,
						gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv()),
						gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())),DABv);
	*/


	//Impongo un residuo e un massimo di iterazioni
	scalar_type oxyres = descr_oxy_transp.Residual_OXY;
	scalar_type err = 1.0; 

	int max_iter = descr_oxy_transp.Max_iterations_OXY;
	int iteration=1;
	
	//salvo la soluzione iniziale in Ct_old e in Cv_old
	gmm::copy(gmm::sub_vector(UM_oxy,
				gmm::sub_interval(0, dof_oxy_transp.Ct())), Ct_old);
	gmm::copy(gmm::sub_vector(UM_oxy,
				gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())), Cv_old);

	//da fare solo nel caso di non linearità
	if(!descr_oxy_transp.HEMOADVECTION && !descr_oxy_transp.REACTION){
		cout<<"Any non-linearities are in the system"<<endl;
	}
	else {

while (RK && iteration<=max_iter && err>oxyres)
{

	if(descr_oxy_transp.REACTION) {
	
	#ifdef M3D1D_VERBOSE_
	cout<<"Assembling Rt in FPM..."<<endl;
	#endif

	gmm::scale(gmm::sub_matrix(AM_oxy,
					gmm::sub_interval(0,dof_oxy_transp.Ct()),
					gmm::sub_interval(0,dof_oxy_transp.Ct())),0.0);

	vector_type ct_guess(dof_oxy_transp.Ct());

	gmm::copy(Ct_old, ct_guess); //Riga di aggiornamento

	//Coefficient for reaction term:	
	vector_type PRESS50(dof_oxy_transp.Ct(), PARAM.real_value("Pm_50")); 
	gmm::scale(PRESS50, param_oxy_transp.alpha_t());

	gmm::add(PRESS50,  ct_guess);

	for (size_type i=0; i<dof_oxy_transp.Ct(); i++){
	ct_guess[i] = 1.0/ct_guess[i];
	}
	gmm::scale(ct_guess, param_oxy_transp.M0());

	vector_type consump_coeff(mf_R.nb_dof());
	getfem::interpolation(mf_oxy_Ct, mf_R, ct_guess, consump_coeff);
	
	//RR
	//Prtché? cosi da avere la dimensione di ct_guess, consump_coeff e press50 identica a quella di UM_oxy(0, Ct) --> utile per il FPM, e 
	//dare come input a asm_tissue_transp il vettore di coefficienti (consump_coeff, da costruire su mf_coeft) con le dimensioni corrette


	getfem::asm_mass_matrix_param(Rt_temp, mimt, mf_oxy_Ct, mf_R, consump_coeff);

	gmm::add(Rt_temp, 
			  gmm::sub_matrix(AM_oxy, 
					gmm::sub_interval(0, dof_oxy_transp.Ct()), 
				 	gmm::sub_interval(0, dof_oxy_transp.Ct())));

	gmm::add(Dt_temp, 
			  gmm::sub_matrix(AM_oxy, 
					gmm::sub_interval(0, dof_oxy_transp.Ct()), 
				 	gmm::sub_interval(0, dof_oxy_transp.Ct())));


	gmm::add(At_temp, 
			  gmm::sub_matrix(AM_oxy, 
					gmm::sub_interval(0, dof_oxy_transp.Ct()), 
				 	gmm::sub_interval(0, dof_oxy_transp.Ct())));

	gmm::add(Btt_temp, 
			  gmm::sub_matrix(AM_oxy, 
					gmm::sub_interval(0, dof_oxy_transp.Ct()), 
				 	gmm::sub_interval(0, dof_oxy_transp.Ct())));

	gmm::clear(Rt_temp);
 	}




 	if(descr_oxy_transp.HEMOADVECTION){

 	#ifdef M3D1D_VERBOSE_
	cout<<"Assembling Ov in FPM..."<<endl;
	#endif

	gmm::scale(gmm::sub_vector(FM_oxy,
					gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())),0.0);
 	vector_type cv_guess (dof_oxy_transp.Cv());

 	gmm::copy(Cv_old, cv_guess); //Riga di aggiornamento

	size_type shift =0;
	size_type shift_h=0;
	
	scalar_type k1;
	k1= param_oxy_transp.N()*param_oxy_transp.MCHC();
	scalar_type k2;
	k2 = pow((param_oxy_transp.Ps_50()*param_oxy_transp.alpha_t()),param_oxy_transp.delta());

	for(size_type i=0; i<nb_branches; ++i)
	{

		vector_type cv_i(mf_Hi[i].nb_dof()); gmm::clear(cv_i);
		vector_type Hi(mf_Hi[i].nb_dof()); gmm::clear(Hi);
		vector_type Uvi(mf_Uvi[i].nb_dof()); gmm::clear(Uvi);
		vector_type psi(mf_Hi[i].nb_dof()); gmm::clear(psi);

	size_type pos=0;
	for (getfem::mr_visitor mrv(mf_oxy_Cv.linked_mesh().region(i)); !mrv.finished(); ++mrv)
		{
		if(pos == 0)
			{
			cv_i[pos] = cv_guess[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0]];
			//cout << " cv_i [" << pos << "] = cv_guess [" << mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0] << "]" <<"="<<cv_guess[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[0]]<< endl;
			pos ++;
			cv_i[pos] = cv_guess[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]];
			//cout << " cv_i [" << pos << "] = cv_guess [" << mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<cv_guess[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;
			}
		else{
			cv_i[pos] = cv_guess[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]];
			//cout << " cv_i [" << pos << "] = cv_guess [" << mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1] << "]" <<"="<<cv_guess[mf_oxy_Cv.ind_basic_dof_of_element(mrv.cv())[1]]<< endl;
			pos ++;
			}
	//ho fatto questo tipo di ciclo perché sul singolo elemento mi prendeva due volte il nodo iniziale: ora prendo entrambi nodi per il primo ciclo e poi per il secondo prendo solo il secondo valore del nodo
		}

		if(i>0) shift_h += mf_Hi[i-1].nb_dof();
		if(i>0) shift += mf_Uvi[i-1].nb_dof();
		
		gmm::add(gmm::sub_vector(UM_HT, 
			gmm::sub_interval(shift_h, mf_Hi[i].nb_dof())), Hi);
		gmm::add(gmm::sub_vector(UM,
			gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof())) ,  Uvi);

		for (size_type j=0; j<mf_Hi[i].nb_dof(); ++j)
		{
			psi[j] = Hi[j]*k1*pow(cv_i[j], param_oxy_transp.delta())/(pow(cv_i[j], param_oxy_transp.delta())+k2); //calcolo dimensionale
			psi[j] = psi[j]/param_oxy_transp.C(); //per adimensionalizzare
		}
	asm_hemoadvection_rhs_network(Ov, mimv, mf_oxy_Cv, mf_coefvi[i], mf_Uvi[i], mf_coefv, mf_Hi[i], Uvi, param.lambdax(i), param.lambday(i), param.lambdaz(i),  param.R(), psi, meshv.region(i));
	}

	gmm::add(Ov, gmm::sub_vector(FM_oxy, 
					gmm::sub_interval(dof_oxy_transp.Ct(),dof_oxy_transp.Cv())));
	} //fine HEMOADVECTION

		//risolvere il nuovo sistema
		assembly_rhs_oxy_transp();

		cout<<"Risolvo il sistema nuovamente"<<endl;
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
}; //esco dal while
} //fine dell'if (presenza o no di non linearità)
	gmm::copy(Ct_old, gmm::sub_vector(UM_oxy, 
					 gmm::sub_interval(0,dof_oxy_transp.Ct())));
	
	gmm::copy(Cv_old, gmm::sub_vector(UM_oxy,
					 gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())));	

	/*time_G=clock()-time_G;
	cout<< "Iterative Process Time = " << ((float)time_G)/CLOCKS_PER_SEC << " s"<< endl;
	SaveResidual.close();*/

	return true;
}; //end of oxygen fixpoint
	

	//Compute the residuals for mass balance at each junction 
	void oxygen_transport3d1d::mass_balance(void){

	#ifdef M3D1D_VERBOSE_		
	cout << " Compute MBD and MBA "   << endl;
	#endif	

	// initialize the MBD and MBA to zero (clear at eac time step)
	for (size_type i=0; i<Jv_oxy_transp.size(); ++i){
		Jv_oxy_transp[i].MBD=0;
		Jv_oxy_transp[i].MBA=0;
	}	

	size_type shift = 0; //counter for branches
	
	for (size_type i=0; i<mf_Uvi.size(); ++i){ // branch loop 	
		if(i>0) shift += mf_Uvi[i-1].nb_dof(); 
		mesh_region &rg_branch = meshv.region(i); // branch region

		for (size_type j=0; j<Jv_oxy_transp.size(); ++j){ //junction loop
			mesh_region &rg_junction = meshv.region(Jv_oxy_transp[j].rg);  // junction region
			// Iterators for all the branches which flow in the junction j
			std::vector<long signed int>::const_iterator bb = Jv_oxy_transp[j].branches.begin();
			std::vector<long signed int>::const_iterator be = Jv_oxy_transp[j].branches.end();
		
			//Check if outflow of branch i is in junction j			
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
	
				//Compute the diffusive flux
				 DIFF = pi* Ri * Ri*(
						UM_oxy[dof_oxy_transp.Ct()+first_C2]-UM_oxy[dof_oxy_transp.Ct()+first_C] )
						/estimate_h(meshv, ii.cv()) ;
	
	
				//Compute the advective fluxes
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


	
	}; // end of mass_balance


 

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

	//Copy solution
	gmm::copy(gmm::sub_vector(UM_oxy, 
		gmm::sub_interval(0, dof_oxy_transp.Ct())),  Ct);
	gmm::copy(gmm::sub_vector(UM_oxy, 
		gmm::sub_interval(dof_oxy_transp.Ct(), dof_oxy_transp.Cv())), Cv);	

	for(size_type a=0;a<dof_oxy_transp.Ct();a++){
		cout<<"Ct["<<a<<"]= "<<Ct[a]<<endl;
	}
	for(size_type b=0;b<dof_oxy_transp.Cv();b++){
		cout<<"Cv["<<b<<"]= "<<Cv[b]<<endl;
	}


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
