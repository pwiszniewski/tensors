#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.h>
#include <hdf5.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

// macro to help check return status of HDF5 functions
#define HDF5_STATUS_CHECK(status) {                 \
  if(status < 0)                                  \
    fprintf(outfile, "%s:%d: Error with HDF5. status code=%d\n", __FILE__, __LINE__, status); \
  }

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace aoints {

extern "C"
int read_options(std::string name, Options &options)
{
    if (name == "AOINTS"|| options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
      options.add_bool("PRINT_INTEGRALS", true);
      /*- Whether to compute two-electron integrals -*/
      options.add_bool("DO_TEI", true);
      // save to a HDF5 file
      options.add_bool("SAVEHDF5", false);
      options.add_str("HDF5_FILENAME", "integrals.h5");
    }

    return true;
}

extern "C"
PsiReturnType aoints(Options &options)
{
    bool print = options.get_bool("PRINT_INTEGRALS");
    int doTei = options.get_bool("DO_TEI");
    bool savehdf5 = options.get_bool("SAVEHDF5");
    std::string filename = options.get_str("HDF5_FILENAME");
    boost::algorithm::to_lower(filename);

    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();

    // Form basis object:
    // Create a basis set parser object.
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    // Construct a new basis set.
    boost::shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, molecule, "BASIS");

    // The integral factory oversees the creation of integral objects
    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory
            (aoBasis, aoBasis, aoBasis, aoBasis));
    
    // N.B. This should be called after the basis has been built, because the geometry has not been
    // fully initialized until this time.
    molecule->print();
    int nbf[] = { aoBasis->nbf() };
    double nucrep = molecule->nuclear_repulsion_energy();

    int nelectrons = 0;
    for(int i=0;i<molecule->natom();i++)
      nelectrons += molecule->true_atomic_number(i);
    
    nelectrons -= molecule->molecular_charge();
    
    std::string sym = molecule->sym_label();
    
    if (sym != "c1"){
	fprintf(outfile, "ERROR: We assume c1 symmetry. Please make it so: add 'symmetry c1' to the molecule specification!");
	return Failure;
    }

    fprintf(outfile, "\n    Nuclear repulsion energy: %16.8f\n\n", nucrep);
    fprintf(outfile, "\n    Number of electrons: %d\n\n", nelectrons);
    fprintf(outfile, "\n    Number of basis functions: %d\n\n", nbf[0]);
    
    // The matrix factory can create matrices of the correct dimensions...
    boost::shared_ptr<MatrixFactory> factory(new MatrixFactory);
    factory->init_with(1, nbf, nbf);

    // Form the one-electron integral objects from the integral factory
    boost::shared_ptr<OneBodyAOInt> sOBI(integral->ao_overlap());
    boost::shared_ptr<OneBodyAOInt> tOBI(integral->ao_kinetic());
    boost::shared_ptr<OneBodyAOInt> vOBI(integral->ao_potential());
    // Form the one-electron integral matrices from the matrix factory
    SharedMatrix sMat(factory->create_matrix("Overlap"));
    SharedMatrix tMat(factory->create_matrix("Kinetic"));
    SharedMatrix vMat(factory->create_matrix("Potential"));
    SharedMatrix hMat(factory->create_matrix("One Electron Ints"));
    // Compute the one electron integrals, telling each object where to store the result
    sOBI->compute(sMat);
    tOBI->compute(tMat);
    vOBI->compute(vMat);

    if(print){
      fprintf(outfile, "\n  One-electron Integrals\n\n");
      sMat->print();
      tMat->print();
      vMat->print();
      // Form h = T + V by first cloning T and then adding V
      hMat->copy(tMat);
      hMat->add(vMat);
      hMat->print();
    }

    // helper variables for HDF5
    hid_t       file_id, group_id, dataset_id, dataspace_id, attribute_id;
    herr_t      status;

    if(savehdf5){
      file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      HDF5_STATUS_CHECK(file_id);
      
      group_id = H5Gcreate(file_id, "/integrals", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      HDF5_STATUS_CHECK(group_id);
      
      dataspace_id = H5Screate(H5S_SCALAR);
      
      attribute_id = H5Acreate (group_id, "nelectrons", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite (attribute_id, H5T_NATIVE_INT, &nelectrons );
      HDF5_STATUS_CHECK(status);
      status = H5Aclose(attribute_id);
      HDF5_STATUS_CHECK(status);
      
      attribute_id = H5Acreate (group_id, "sp_dim", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite (attribute_id, H5T_NATIVE_INT, &nbf[0] );
      HDF5_STATUS_CHECK(status);
      status = H5Aclose(attribute_id);
      HDF5_STATUS_CHECK(status);
      
      attribute_id = H5Acreate (group_id, "nuclear_repulsion_energy", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite (attribute_id, H5T_NATIVE_DOUBLE, &nucrep );
      HDF5_STATUS_CHECK(status);
      status = H5Aclose(attribute_id);
      HDF5_STATUS_CHECK(status);

      
      status = H5Sclose(dataspace_id);
      HDF5_STATUS_CHECK(status);

      // Save the size of irrep as 1-1 array (Octave HDF5 parser is ugly)
      hsize_t dimarray = 1;
      dataspace_id = H5Screate_simple(1, &dimarray, NULL);

      dataset_id = H5Dcreate(group_id, "Nbf", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nbf[0] );
      HDF5_STATUS_CHECK(status);
      status = H5Dclose(dataset_id);
      HDF5_STATUS_CHECK(status);
      
      status = H5Sclose(dataspace_id);
      HDF5_STATUS_CHECK(status);
      /*
      hsize_t dimarray = nbf[0] * nbf[0];
      dataspace_id = H5Screate_simple(1, &dimarray, NULL);

      dataset_id = H5Dcreate(group_id, "Overlap", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, sMat->get_pointer(0) );
      HDF5_STATUS_CHECK(status);
      status = H5Dclose(dataset_id);
      HDF5_STATUS_CHECK(status);
      
      dataset_id = H5Dcreate(group_id, "Kinetic", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tMat->get_pointer(0) );
      HDF5_STATUS_CHECK(status);
      status = H5Dclose(dataset_id);
      HDF5_STATUS_CHECK(status);
      
      dataset_id = H5Dcreate(group_id, "Potential", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vMat->get_pointer(0) );
      HDF5_STATUS_CHECK(status);
      status = H5Dclose(dataset_id);
      HDF5_STATUS_CHECK(status);
      
      
      status = H5Sclose(dataspace_id);
      HDF5_STATUS_CHECK(status);
      */
    }
    
    if(doTei){
      fprintf(outfile, "\n  Two-electron Integrals\n\n");
      
      shared_ptr<Matrix> VMat;
      
      if(savehdf5){
	hsize_t dimarray = nbf[0] * nbf[0] * nbf[0] * nbf[0];
	dataspace_id = H5Screate_simple(1, &dimarray, NULL);
	
	VMat.reset(new Matrix(nbf[0] * nbf[0], nbf[0] * nbf[0]));
	VMat->set(0.0);
      }
      
      // Now, the two-electron integrals
      boost::shared_ptr<TwoBodyAOInt> eri(integral->eri());
      // The buffer will hold the integrals for each shell, as they're computed
      const double *buffer = eri->buffer();
      // The iterator conveniently lets us iterate over functions within shells
      AOShellCombinationsIterator shellIter = integral->shells_iterator();
      int count=0;
      for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
	// Compute quartet
	eri->compute_shell(shellIter);
            // From the quartet get all the integrals
	AOIntegralsIterator intIter = shellIter.integrals_iterator();
	for (intIter.first(); intIter.is_done() == false; intIter.next()) {
	  int p = intIter.i();
	  int q = intIter.j();
	  int r = intIter.k();
	  int s = intIter.l();
	  
	  //	  if(print)
	  //  fprintf(outfile, "\t(%2d %2d | %2d %2d) = %20.15f\n",
	  //	    p, q, r, s, buffer[intIter.index()]);

	  // PSI is chemical notation, we get: (pr|V|qs)
	  // so: V_abcd = V(a*nmo+b,c*nmo+d)

	  if(savehdf5){
	    //Cumulant order
	    int idx1 = p * nbf[0] + q;
	    int idx2 = r * nbf[0] + s;
	    VMat->set(idx1, idx2, buffer[intIter.index()]);
	    VMat->set(idx2, idx1, buffer[intIter.index()]);
	    
	    idx1 = p * nbf[0] + r;
	    idx2 = q * nbf[0] + s;
	    VMat->set(idx1, idx2, buffer[intIter.index()]);
	    VMat->set(idx2, idx1, buffer[intIter.index()]);
	    
	    idx1 = s * nbf[0] + r;
	    idx2 = q * nbf[0] + p;
	    VMat->set(idx1, idx2, buffer[intIter.index()]);
	    VMat->set(idx2, idx1, buffer[intIter.index()]);

	    idx1 = s * nbf[0] + q;
	    idx2 = r * nbf[0] + p;
	    VMat->set(idx1, idx2, buffer[intIter.index()]);
	    VMat->set(idx2, idx1, buffer[intIter.index()]);

	    /*
	    int idx1 = p * nbf[0] + r;
	    int idx2 = q * nbf[0] + s;
	    VMat->set(idx1, idx2, buffer[intIter.index()]);
	    VMat->set(idx2, idx1, buffer[intIter.index()]);
	    
	    idx1 = r * nbf[0] + p;
	    idx2 = s * nbf[0] + q;
	    VMat->set(idx1, idx2, buffer[intIter.index()]);
	    VMat->set(idx2, idx1, buffer[intIter.index()]);
	    
	    // with real orbitals we can swap r and s
	    idx1 = q * nbf[0] + r;
	    idx2 = p * nbf[0] + s;
	    VMat->set(idx1, idx2, buffer[intIter.index()]);
	    VMat->set(idx2, idx1, buffer[intIter.index()]);

	    idx1 = r * nbf[0] + q;
	    idx2 = s * nbf[0] + p;
	    VMat->set(idx1, idx2, buffer[intIter.index()]);
	    VMat->set(idx2, idx1, buffer[intIter.index()]);
	    */
	    //Mulliken order
	    /* int idx1 = p * nbf[0] + q;
	    int idx2 = r * nbf[0] + s;
	    VMat->set(idx1, idx2, buffer[intIter.index()]);
	    VMat->set(idx2, idx1, buffer[intIter.index()]);
	    
	    idx1 = q * nbf[0] + p;
	    idx2 = s * nbf[0] + r;
	    VMat->set(idx1, idx2, buffer[intIter.index()]);
	    VMat->set(idx2, idx1, buffer[intIter.index()]);
	    
	    // with real orbitals we can swap r and s
	    idx1 = r * nbf[0] + q;
	    idx2 = p * nbf[0] + s;
	    VMat->set(idx1, idx2, buffer[intIter.index()]);
	    VMat->set(idx2, idx1, buffer[intIter.index()]);

	    idx1 = q * nbf[0] + r;
	    idx2 = s * nbf[0] + p;
	    VMat->set(idx1, idx2, buffer[intIter.index()]);
	    VMat->set(idx2, idx1, buffer[intIter.index()]);
	    */
	  }

	  ++count;
	}
      }
      fprintf(outfile, "\n\tThere are %d unique integrals\n\n", count);

      if(savehdf5){
	dataset_id = H5Dcreate(group_id, "TwoBody", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, VMat->get_pointer(0) );
	HDF5_STATUS_CHECK(status);
	status = H5Dclose(dataset_id);
	HDF5_STATUS_CHECK(status);
	
	status = H5Sclose(dataspace_id);
	HDF5_STATUS_CHECK(status);
      }
    }

    if(savehdf5){
      status = H5Gclose(group_id);
      HDF5_STATUS_CHECK(status);
      
      status = H5Fclose(file_id);
      HDF5_STATUS_CHECK(status);
    }
    
    return Success;
}
    
}} // End Namespaces
