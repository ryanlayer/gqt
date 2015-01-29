#ifndef __PLT_H__
#define __PLT_H__

#include "genotq.h"

struct plt_file {
    FILE *file;
    uint32_t num_fields, num_records;
    long header_offset;
};

struct plt_file init_plt_file(char *file_name);

/**
 * @brief Convert a VCF file to a by-variant plain text file
 *
 * @param in_file_name VCF file name
 * @param num_fields Number of fields in the VCF
 * @param num_records Number of records in the VCF
 * @param out_file_name plain text file name
 *
 * @retval 0 if things worked
 * @retval 1 otherwise
 *
 * Example Usage:
 * @code
 * @endcode
 */
int convert_file_by_name_vcf_to_plt(char *in_file_name,
                                    uint32_t num_fields,
                                    uint32_t num_records,
                                    char *out_file_name);

/**
 * @brief Switch a plain text records to fields
 *
 * @param in_file_name plain text input file name
 * @param out_file_name plain text output file name
 *
 * @retval 0 if things worked
 * @retval 1 otherwise
 *
 * Example Usage:
 * @code
 * @endcode
 */
int convert_file_by_name_invert_plt(char *in_file_name, char *out_file_name);

/**
 * @brief Convert a plain text by-variant file to a VCF file
 *
 * @param in_file_name Plain text file name
 * @param out_file_name VCF file name
 *
 * @retval 1 if things worked
 * @retval 0 otherwise
 *
 * Example Usage:
 * @code
 *     char *plt_file_name="data/10.1e4.var.txt";
 *     char *vcf_file_name="data/10.1e4.var.vcf";
 *
 *     int r = convert_file_by_name_plt_to_vcf(plt_file_name,vcf_file_name);
 * @endcode
 */
int convert_file_by_name_plt_to_vcf(char *in_file_name, char *out_file_name);

/**
 * @brief Convert a plain text by-variant file to a VCF file
 *
 * @param pf initialized by-variant plt file data structure
 * @param out_file_name name of the VCF file
 *
 * @retval 1 if things worked
 * @retval 0 otherwise
 *
 * Example Usage:
 * @code
 *     char *plt_file_name="data/10.1e4.var.txt";
 *     char *vcf_file_name="data/10.1e4.var.vcfn";
 *
 *     struct plt_file pf = init_plt_file(plt_file_name);
 *     int r = convert_file_plt_to_vcf(pf, vcf_file_name);
 *     fclose(pf.file);
 * @endcode
 */
int convert_file_plt_to_vcf(struct plt_file pf, char *out_file_name);


/**
 * @brief Convert a plain text file to an uncompressed binary file
 *
 * @param in_file_name Plain text file name
 * @param out_file_name Uncompressed binary file name
 *
 * @retval 1 if things worked
 * @retval 0 otherwise
 *
 * Example Usage:
 * @code
 *     char *plt_file_name="data/10.1e4.ind.txt";
 *     char *ubin_file_name="data/10.1e4.ind.ubin";
 *
 *     int r = convert_file_by_name_plt_to_ubin(plt_file_name,ubin_file_name);
 * @endcode
 */
int convert_file_by_name_plt_to_ubin(char *in_file_name, char *out_file_name);

/**
 * @brief Convert a plain text file to an uncompressed binary file
 *
 * @param pf initialized plt file data structure
 * @param out_file_name name of the uncompressed binary file
 *
 * @retval 1 if things worked
 * @retval 0 otherwise
 *
 * Example Usage:
 * @code
 *     char *plt_file_name="data/10.1e4.ind.txt";
 *     char *ubin_file_name="data/10.1e4.ind.ubin";
 *
 *     struct plt_file pf = init_plt_file(plt_file_name);
 *     int r = convert_file_plt_to_ubin(pf, ubin_file_name);
 *     fclose(pf.file);
 * @endcode
 */
int convert_file_plt_to_ubin(struct plt_file pf, char *out_file_name);

uint32_t plt_line_to_packed_ints(char *line,
                                     uint32_t num_fields, 
                                     uint32_t **packed_ints);
/**
 * @brief Convert an string of plain text values over the alphabet [0,1,2,3]
 *        to WAH encoding of the bitmap index representation of the two-bit
 *        binary
 *
 *  In this scheme, there are 4 uniq values in plain text (0,1,2,3), so the WAH
 *  encoding will create 4 different sets, one for each bitmap index
 *
 * @param plt         Space-seperated string of 0,1,2, or 3
 * @param plt_len     Number of digits in P
 * @param W         The resulting array of WAH words (memory allocted in 
 *                  function, value set in function)
 * @param wah_sizes   A list the 4 sizes
 *
 * @retval number of ints in W
 *
 * Example Usage:
 * @code
 *      char *plt = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
 *                  "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
 *                  "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
 *                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
 *                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
 *                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
 *                  "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
 *                  "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
 *                  "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
 *                  "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 "
 *                  "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 "
 *                  "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3";
 *
 *      uint32_t *wah;
 *      uint32_t *wah_sizes;
 *      uint32_t wah_len = plt_to_bitmap_wah(plt,
 *                                               192,
 *                                               &wah,
 *                                               &wah_sizes);
 * @endcode
 */
uint32_t plt_to_bitmap_wah(char *plt,
                               uint32_t plt_len,
                               uint32_t **W,
                               uint32_t **wah_sizes);

/**
 * @brief Get a pointer to the packed int representation of a particular plain
 * text record
 *
 * @param pf The initialized plain text file data structure
 * @param plt_record The record ID
 * @param plt A pointer set within the fuction that points to the record
 *            of intrest
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 *      char *plt_file_name="data/10.1e4.ind.txt";
 *      struct plt_file pf = init_wah_file(plt_file_name);
 *      uint32_t *plt;
 *      unsinged int plt_size = get_plt_record(pf,
 *                                             test_record,
 *                                             &plt);
 * @endcode
 */
uint32_t get_plt_record(struct plt_file pf,
                            uint32_t plt_record,
                            uint32_t **plt);

/**
 * @brief Return a count of records whose values are >= start_test_value and <
 * end_test_value
 *
 * @param pf The initialized plain text encoded file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R array of ints that contain counts
 *
 * @retval number of ints in R
 *
 * Example Usage:
 * @code
 * @endcode
 */
 uint32_t count_range_records_plt(struct plt_file pt,
                                      uint32_t *record_ids,
                                      uint32_t num_r,
                                      uint32_t start_test_value,
                                      uint32_t end_test_value,
                                      uint32_t **R);
/**
 * @brief Return records whose values are >= start_test_value and <
 * end_test_value
 *
 * @param pf The initialized plain text encoded file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R packed ints where the bit is zero if all fields meet that condition
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
 uint32_t range_records_plt(struct plt_file pt,
                              uint32_t *record_ids,
                              uint32_t num_r,
                              uint32_t start_test_value,
                              uint32_t end_test_value,
                              uint32_t **R);
/**
 * @brief Return fields whose values are >= start_test_value and <
 * end_test_value
 *
 * @param pf The initialized plain text encoded file
 * @param field_ids array of integer ids of the fields to test
 * @param num_f number of fields in fieldids
 * @param test_value value to test fields against
 * @param R packed ints where the bit is zero if all fields meet that condition
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
 uint32_t range_fields_plt(struct plt_file pt,
                               uint32_t *field_ids,
                               uint32_t num_f,
                               uint32_t start_test_value,
                               uint32_t end_test_value,
                               uint32_t **R);

/*
 * @brief Return records whose value is equal to the test value
 *
 * @param pf The initialized plain text encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
uint32_t eq_records_plt(struct plt_file pf,
                            uint32_t *record_ids,
                            uint32_t num_r,
                            uint32_t test_value,
                            uint32_t **R);

/**
 * @brief Return records whose value is NOT equal to the test value
 *
 * @param pf The initialized plain text encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
uint32_t ne_records_plt(struct plt_file pf,
                            uint32_t *record_ids,
                            uint32_t num_r,
                            uint32_t test_value,
                            uint32_t **R);

/**
 * @brief Return records whose value are greater than the test value
 *
 * @param pf The initialized plain text encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
uint32_t gt_records_plt(struct plt_file pf,
                            uint32_t *record_ids,
                            uint32_t num_r,
                            uint32_t test_value,
                            uint32_t **R);

/**
 * @brief Return fields whose value are greater than the test value
 *
 * @param pf The initialized plain text encoded bitmap file
 * @param fields_ids array of integer ids of the fields to test
 * @param num_f number of fields in field_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
uint32_t gt_fields_plt(struct plt_file pf,
                           uint32_t *field_ids,
                           uint32_t num_f,
                           uint32_t test_value,
                           uint32_t **R);


/**
 * @brief Return filds whose value are greater than the test value
 *
 * @param pf The initialized plain text encoded bitmap file
 * @param field_ids array of integer ids of the fields to test
 * @param num_f number of fields in field_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
uint32_t gt_records_plt(struct plt_file pf,
                            uint32_t *field_ids,
                            uint32_t num_f,
                            uint32_t test_value,
                            uint32_t **R);


/**
 * @brief Return the count of records whose value are greater than the test
 * value
 *
 * @param pf The initialized plain text encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R array of integers with counts
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
uint32_t gt_records_plt(struct plt_file pf,
                            uint32_t *record_ids,
                            uint32_t num_r,
                            uint32_t test_value,
                            uint32_t **R);



/**
 * @brief Return records whose value are greater than or equal to the test value
 *
 * @param pf The initialized plain text encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
uint32_t gte_records_plt(struct plt_file pf,
                            uint32_t *record_ids,
                            uint32_t num_r,
                            uint32_t test_value,
                            uint32_t **R);


/**
 * @brief Return records whose value are less than the test value
 *
 * @param pf The initialized plain text encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
uint32_t lt_records_plt(struct plt_file pf,
                            uint32_t *record_ids,
                            uint32_t num_r,
                            uint32_t test_value,
                            uint32_t **R);

/**
 * @brief Return records whose value are less or equal to the test value
 *
 * @param pf The initialized plain text encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */
uint32_t lte_records_plt(struct plt_file pf,
                            uint32_t *record_ids,
                            uint32_t num_r,
                            uint32_t test_value,
                            uint32_t **R);

/**
 * @brief Return the number of records whose value are greater than the test
 * value
 *
 * @param pf The initialized plain text encoded bitmap file
 * @param record_ids array of integer ids of the records to test
 * @param num_r number of records in record_ids
 * @param test_value value to test fields against
 * @param R result with 
 *
 * @retval number of ints in the record
 *
 * Example Usage:
 * @code
 * @endcode
 */

uint32_t gt_count_records_plt(struct plt_file pf,
                                  uint32_t *record_ids,
                                  uint32_t num_r,
                                  uint32_t test_value,
                                  uint32_t **R);

/**
 * @brief Print a plaint text file.
 *
 * If num_r > 0, then record_ids should contain the ids of records in the file
 * that will be displayed, otherwise all records will be displayed.
 *
 * @param pf An intitialized plaint text file
 * @param record_ids An array of record ids
 * @param num_r number of records in the array
 *
 * @retval number of records printed 
 */
uint32_t print_plt(struct plt_file pf,
                       uint32_t *record_ids,
                       uint32_t num_r);

/**
 * @brief Print plain text file (wrapper around print_plt). 
 *
 * If num_r > 0, then record_ids should contain the ids of records in the file
 * that will be displayed, otherwise all records will be displayed.
 *
 * @param pf_file_name Plain text file name
 * @param record_ids An array of record ids
 * @param num_r number of records in the array
 *
 * @retval number of records printed 
 */
uint32_t print_by_name_plt(char *pf_file_name,
                               uint32_t *record_ids,
                               uint32_t num_r);

/**
 * @brief 
 *
 * @param ints
 * @param num_ints
 *
 * @retval 
 */

uint32_t pack_2_bit_ints(int *ints, int num_ints);

/**
 * @brief Invert a plain text line into an array of uncompressed binary encoded
 * values
 *
 * The function is used to build the inverted uncompressed binary values one
 * plain text line at a time.  The ubin value must be set to NULL to start with
 * so that this function can allocate the neccessary space.
 *
 * @param line              A plain text line of genotypes
 * @param num_fields        Number of fields in the line
 * @param num_records       Total number of records
 * @param new_num_fields    The number of fields in the uncompressed binary
 * @param new_num_records   The number of records in the uncompressed binary
 * @param field_i           The current field index within the uncompressed
 * binary 
 * @param ubin              The uncompressed binary values
 *
 * @retval The next field index in the uncompressed binary 
 */
uint32_t invert_plt_to_ubin(char *line,
                                uint32_t num_fields,
                                uint32_t num_records,
                                uint32_t *new_num_fields,
                                uint32_t *new_num_records,
                                uint32_t field_i,
                                uint32_t ***ubin);
/**
 * @brief Invert a plain text file and store as an uncompressed binary
 *
 * @param in_file_name Input plain text file
 * @param out_file_name Output uncompressed binary file
 *
 * @retval 0 if everything worked
 */
int convert_file_by_name_invert_plt_to_ubin(char *in_file_name,
                                            char *out_file_name);
#endif
