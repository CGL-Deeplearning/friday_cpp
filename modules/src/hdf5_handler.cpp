//
// Created by Kishwar Shafin on 6/25/18.
//
#include "../headers/hdf5_handler.h"

hdf5_handler::hdf5_handler(const string file_name) {
    this->file_name = file_name;
    hid_t   file_id;
    file_id = H5Fcreate (file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5Fclose (file_id);
}

void hdf5_handler::save_img_to_hdf5(const int candidate_image[IMAGE_HEIGHT][IMAGE_WIDTH][TOTAL_CHANNELS]) {
    hid_t       dataset_id, dataspace_id;  /* identifiers */
    herr_t      status;
    hid_t       plist_id;
    hid_t       file_id;
    file_id = H5Fopen(this->file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    /* Create dataset "Compressed Data" in the group using absolute name.  */
    const hsize_t dims[] = { IMAGE_HEIGHT, IMAGE_WIDTH, TOTAL_CHANNELS };
    dataspace_id = H5Screate_simple (3, dims, NULL);

    plist_id  = H5Pcreate (H5P_DATASET_CREATE);

    /* Dataset must be chunked for compression */
    const hsize_t cdims[] = { 20, 20, 5};
    status = H5Pset_chunk (plist_id, 3, cdims);

    status = H5Pset_deflate (plist_id, 6);

    hid_t mem_space = H5Screate_simple(3, dims, NULL);

    dataset_id = H5Dcreate2 (file_id, "image", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);

    status = H5Dwrite (dataset_id, H5T_NATIVE_INT, mem_space, H5S_ALL, H5P_DEFAULT, candidate_image);

    status = H5Sclose (dataspace_id);
    status = H5Dclose (dataset_id);
    status = H5Pclose (plist_id);
    status = H5Fclose (file_id);
}


void hdf5_handler::save_allele_info(string chromosome_name, string position,
                                    string allele1, string allele1_type,
                                    string allele2, string allele2_type) {

    hid_t   file_id;
    file_id = H5Fopen(this->file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    const char  *wdata[6] = {chromosome_name.c_str(), position.c_str(),
                             allele1.c_str(), allele1_type.c_str(),
                             allele2.c_str(), allele2_type.c_str()};

    hid_t       filetype, memtype, space, dset;
    herr_t      status;
    hsize_t     dims[1] = {6};

    filetype = H5Tcopy (H5T_FORTRAN_S1);
    status = H5Tset_size (filetype, H5T_VARIABLE);
    memtype = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (memtype, H5T_VARIABLE);

    /*
     * Create dataspace.  Setting maximum size to NULL sets the maximum
     * size to be the current size.
     */
    space = H5Screate_simple (1, dims, NULL);
    /*
     * Create the dataset and write the variable-length string data to
     * it.
     */
    dset = H5Dcreate (file_id, "allele_info", filetype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
    /*
     * Close and release resources.
     */
    status = H5Dclose (dset);
    status = H5Sclose (space);
    status = H5Tclose (filetype);
    status = H5Tclose (memtype);
    status = H5Fclose (file_id);

}

hdf5_handler::~hdf5_handler() {

}