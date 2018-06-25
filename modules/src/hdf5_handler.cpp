//
// Created by Kishwar Shafin on 6/25/18.
//
#include "../headers/hdf5_handler.h"

void hdf5_handler::save_img_to_hdf5(image_generator candidate_image) {
    hid_t       file_id, dataset_id, dataspace_id;  /* identifiers */
    herr_t      status;
    hid_t       plist_id;

    file_id = H5Fcreate ("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

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

    status = H5Dwrite (dataset_id, H5T_NATIVE_INT, mem_space, H5S_ALL, H5P_DEFAULT, candidate_image.image_array);

    status = H5Sclose (dataspace_id);
    status = H5Dclose (dataset_id);
    status = H5Pclose (plist_id);
    status = H5Fclose (file_id);
}