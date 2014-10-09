#include "genotq.h"
#include <hdf5.h>
#include <assert.h>
#include <string.h>

//{{{struct hdf5_file init_hdf5_file(char *hdf5_file_name,
struct hdf5_file init_hdf5_file(char *hdf5_file_name,
                                uint32_t num_vars)
{

    struct hdf5_file hdf5_f;
    hdf5_f.num_gt_ints = 1 + ((num_vars - 1) / 32);

    hdf5_f.file_id = H5Fcreate(hdf5_file_name,
                               H5F_ACC_TRUNC,
                               H5P_DEFAULT,
                               H5P_DEFAULT);

    assert(hdf5_f.file_id >= 0);

    hdf5_f.plist_id = H5Pcreate(H5P_DATASET_CREATE);

    assert(hdf5_f.plist_id >= 0);

    hsize_t  c_dim[1];
    if (hdf5_f.num_gt_ints < 10)
        c_dim[0] = hdf5_f.num_gt_ints;
    else
        c_dim[0] = 10;

    int status = H5Pset_chunk(hdf5_f.plist_id, 1, c_dim);
    assert(status >= 0);

    status = H5Pset_deflate (hdf5_f.plist_id, 6);
    assert(status >= 0);

    hsize_t gt_dim[1];
    gt_dim[0] = hdf5_f.num_gt_ints;

    hdf5_f.gt_dataspace_id = H5Screate_simple(1, gt_dim, NULL);
    assert(hdf5_f.gt_dataspace_id >= 0);


    hsize_t md_dim[1];
    md_dim[0] = 1;
    hdf5_f.md_dataspace_id = H5Screate_simple(1, md_dim, NULL);
    assert(hdf5_f.md_dataspace_id >= 0);

    return hdf5_f;
}
//}}}

//{{{int write_hdf5_gt(struct hdf5_file hdf5_f,
int write_hdf5_gt(struct hdf5_file hdf5_f,
                  uint32_t id,
                  uint32_t *gt,
                  char *md)
{
    /* write the genotype data */
    char gt_name[10];
    sprintf(gt_name, "%u", id);
    hid_t dataset_id = H5Dcreate2(hdf5_f.file_id,
                                  gt_name,
                                  H5T_STD_I32BE,
                                  hdf5_f.gt_dataspace_id,
                                  H5P_DEFAULT,
                                  hdf5_f.plist_id,
                                  H5P_DEFAULT);
    assert(dataset_id >= 0);

    int status = H5Dwrite(dataset_id,
                          H5T_NATIVE_INT,
                          H5S_ALL,
                          H5S_ALL,
                          H5P_DEFAULT,
                          gt);
    assert(status >= 0);

    status = H5Dclose(dataset_id);
    assert(status >= 0);

    /* write the variant metadata data */

    /* Create a data type to refer to. */
    hid_t datatype = H5Tcopy (H5T_C_S1);
    assert(datatype >= 0);

    status = H5Tset_size (datatype, H5T_VARIABLE);
    assert(status >= 0);

    /* create the datatype properties */
    hid_t props = H5Pcreate (H5P_DATASET_CREATE);
    assert (props >= 0);

    char md_name[11];
    sprintf(md_name, "m%u", id);
    dataset_id = H5Dcreate(hdf5_f.file_id,
                           md_name,
                           datatype,
                           hdf5_f.md_dataspace_id,
                           H5P_DEFAULT,
                           props,
                           H5P_DEFAULT);
    assert(dataset_id >= 0);

    char *md_a[1];
    md_a[0] = md;

    status = H5Dwrite(dataset_id,
                      datatype,
                      H5S_ALL,
                      H5S_ALL,
                      H5P_DEFAULT,
                      md_a);
    assert(status >= 0);

    status = H5Dclose(dataset_id);
    assert(status >= 0);

    return status;
}
//}}}

//{{{int read_hdf5_gt(struct hdf5_file hdf5_f, uint32_t id, uint32_t *gt)
int read_hdf5_gt(struct hdf5_file hdf5_f, uint32_t id, uint32_t *gt)
{
    char name[10];
    sprintf(name, "%u", id);
    hid_t dataset_id = H5Dopen2(hdf5_f.file_id,
                                name,
                                H5P_DEFAULT);

    int status = H5Dread(dataset_id,
                         H5T_NATIVE_INT,
                         H5S_ALL,
                         H5S_ALL,
                         H5P_DEFAULT,
                         gt);
    status = H5Dclose (dataset_id);
    return status;
}
//}}}

//{{{int read_hdf5_md(struct hdf5_file hdf5_f, uint32_t id, char **md)
int read_hdf5_md(struct hdf5_file hdf5_f, uint32_t id, char **md)
{
    char md_name[11];
    sprintf(md_name, "m%u", id);
    hid_t dataset_id = H5Dopen2(hdf5_f.file_id,
                                md_name,
                                H5P_DEFAULT);


    hid_t datatype = H5Dget_type(dataset_id);
    /* Create a data type to refer to. */
    //hid_t datatype = H5Tcopy (H5T_C_S1);
    
    char *out[0];

    int status = H5Dread(dataset_id,
                         datatype,
                         H5S_ALL,
                         H5S_ALL,
                         H5P_DEFAULT,
                         out);

    *md = (char *) malloc(strlen(out[0]) * sizeof(char));
    strcpy(*md, out[0]);
    free(out[0]);

    status = H5Dclose (dataset_id);
    status = H5Tclose (datatype);

    
    return status;
}
//}}}

//{{{int close_hdf5_file(struct hdf5_file hdf5_f)
int close_hdf5_file(struct hdf5_file hdf5_f)
{
    int status = H5Sclose(hdf5_f.gt_dataspace_id);
    status |= H5Pclose(hdf5_f.plist_id);
    status |= H5Fclose(hdf5_f.file_id);

    return status;
}
//}}}
