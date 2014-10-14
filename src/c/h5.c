#include "genotq.h"
#include <stdio.h>
#include <hdf5.h>
#include <assert.h>
#include <string.h>
#include <sys/param.h>

//{{{struct hdf5_file init_hdf5_file(char *hdf5_file_name,
struct hdf5_file init_hdf5_file(char *hdf5_file_name,
                                uint32_t num_vars,
                                uint32_t num_inds)
{
    struct hdf5_file hdf5_f;

    hdf5_f.num_vars = num_vars;
    hdf5_f.num_inds = num_inds;

    hdf5_f.num_gt_ints = 1 + ((hdf5_f.num_inds - 1) / 16);
    hdf5_f.num_r_gt_ints = 1 + ((hdf5_f.num_vars - 1) / 16);

    // Create the file
    hdf5_f.file_id = H5Fcreate(hdf5_file_name,
                               H5F_ACC_TRUNC,
                               H5P_DEFAULT,
                               H5P_DEFAULT);

    assert(hdf5_f.file_id >= 0);

    // Set up compression
    hdf5_f.plist_id = H5Pcreate(H5P_DATASET_CREATE);
    assert(hdf5_f.plist_id >= 0);

    hsize_t  c_dim[1];
    if (MIN(hdf5_f.num_gt_ints, hdf5_f.num_r_gt_ints)  < 10)
        c_dim[0] = MIN(hdf5_f.num_gt_ints, hdf5_f.num_r_gt_ints);
    else
        c_dim[0] = 10;

    int status = H5Pset_chunk(hdf5_f.plist_id, 1, c_dim);
    assert(status >= 0);

    status = H5Pset_deflate (hdf5_f.plist_id, 6);
    assert(status >= 0);

    // Set up genotype space
    hsize_t gt_dim[1];
    gt_dim[0] = hdf5_f.num_gt_ints;

    hdf5_f.gt_dataspace_id = H5Screate_simple(1, gt_dim, NULL);
    assert(hdf5_f.gt_dataspace_id >= 0);


    // Set up metadata space
    hsize_t md_dim[1];
    md_dim[0] = 1;
    hdf5_f.md_dataspace_id = H5Screate_simple(1, md_dim, NULL);
    assert(hdf5_f.md_dataspace_id >= 0);

    // Set up rotated genotype space
    hsize_t r_gt_dim[1];
    r_gt_dim[0] = hdf5_f.num_r_gt_ints;

    hdf5_f.r_gt_dataspace_id = H5Screate_simple(1, r_gt_dim, NULL);
    assert(hdf5_f.r_gt_dataspace_id >= 0);


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

//{{{int read_hdf5_r_gt(struct hdf5_file hdf5_f, uint32_t id, uint32_t *gt)
int read_hdf5_r_gt(struct hdf5_file hdf5_f, uint32_t id, uint32_t *r_gt)
{
    char name[10];
    sprintf(name, "r%u", id);
    hid_t dataset_id = H5Dopen2(hdf5_f.file_id,
                                name,
                                H5P_DEFAULT);

    int status = H5Dread(dataset_id,
                         H5T_NATIVE_INT,
                         H5S_ALL,
                         H5S_ALL,
                         H5P_DEFAULT,
                         r_gt);
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

//{{{int 
int init_r_gt(struct hdf5_file hdf5_f)
{
    uint32_t *r_gt = (uint32_t *) 
        malloc(hdf5_f.num_r_gt_ints * sizeof(uint32_t));

    uint32_t i;
    char r_gt_name[11];
    for (i = 0; i < hdf5_f.num_inds; ++i) {
        sprintf(r_gt_name, "r%u", i);
        hid_t dataset_id = H5Dcreate2(hdf5_f.file_id,
                                      r_gt_name,
                                      H5T_STD_I32BE,
                                      hdf5_f.r_gt_dataspace_id,
                                      H5P_DEFAULT,
                                      hdf5_f.plist_id,
                                      H5P_DEFAULT);
        assert(dataset_id >= 0);

        int status = H5Dwrite(dataset_id,
                              H5T_NATIVE_INT,
                              H5S_ALL,
                              H5S_ALL,
                              H5P_DEFAULT,
                              r_gt);
        assert(status >= 0);

        status = H5Dclose(dataset_id);
        assert(status >= 0);
    }

    free(r_gt);

    return 0;

}
//}}}

//{{{int set_r_gt(struct hdf5_file hdf5_f,
int set_r_gt(struct hdf5_file hdf5_f,
             uint32_t r_gt_i,
             uint32_t i,
             uint32_t v)
{
    hsize_t offset[1] = {i};
    hsize_t count[1] = {1};

    char r_gt_name[11];
    sprintf(r_gt_name, "r%u", r_gt_i);

    hid_t dataset = H5Dopen(hdf5_f.file_id, r_gt_name, H5P_DEFAULT);
    hid_t dataspace = H5Dget_space (dataset); 


    int status = H5Sselect_hyperslab (dataspace,
                                      H5S_SELECT_SET,
                                      offset,
                                      NULL, 
                                      count,
                                      NULL);

    hsize_t dim[1] = {1};  
    hid_t memspace = H5Screate_simple (1, dim, NULL);   

    uint32_t o[1] = {v};

    status = H5Dwrite (dataset,
                       H5T_NATIVE_INT,
                       memspace,
                       dataspace,
                       H5P_DEFAULT,
                       o);

    return status;
}
//}}}

void sort_rotate_gt_md(pri_queue *q,
                       struct hdf5_file *hdf5_f,
                       char *md_out_file_name)
{

    FILE *f = fopen(md_out_file_name, "w");

    uint32_t i, two_bit_i = 0, int_i = 0;
    int r = init_r_gt(*hdf5_f);

    char *md_out;

    uint32_t *gt_out =
            (uint32_t *)malloc(hdf5_f->num_gt_ints * sizeof(uint32_t));

    uint32_t *gt_buff =
            (uint32_t *)malloc(hdf5_f->num_inds * sizeof(uint32_t));
    memset(gt_buff, 0 , hdf5_f->num_inds * sizeof(uint32_t));

    priority p;

    while ( priq_top(*q, &p) != NULL ) {
        int *d = priq_pop(*q, &p);

        read_hdf5_md(*hdf5_f, *d, &md_out);
        fprintf(f, "%s\n", md_out);
        free(md_out);

        read_hdf5_gt(*hdf5_f, *d, gt_out);

        for (i = 0; i < hdf5_f->num_inds; ++i) {
            //printf("%u ", (gt_out[0] >> (30 - i*2)) & 3);
            int bb = (gt_out[0] >> (30 - i*2)) & 3;
            gt_buff[i] += bb << (30 - two_bit_i*2);
        }
        //printf("\n");

        two_bit_i += 1;
        if (two_bit_i == 16) {
            for (i = 0; i < hdf5_f->num_inds; ++i) {
                set_r_gt(*hdf5_f, i, int_i, gt_buff[i]);
            }

            memset(gt_buff, 0, hdf5_f->num_inds * sizeof(uint32_t));
            two_bit_i = 0;
            int_i += 1;
        }

    }

    if (two_bit_i > 0) {
        for (i = 0; i < hdf5_f->num_inds; ++i) {
            set_r_gt(*hdf5_f, i, int_i, gt_buff[i]);
        }
    }

    free(gt_out);
    free(gt_buff);
    fclose(f);
}



