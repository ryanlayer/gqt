#ifndef __UBIN_H__
#define __UBIN_H__

#include <stdint.h>

struct ubin_file {
    FILE *file;
    char *file_name;
    uint32_t num_fields, num_records;
    long header_offset;
};

struct ubin_file init_ubin_file(char *file_name);

/**
 * @brief Convert an array of uncompressed binary values to a bitmap index of
 *        the values.
 *
 * @param U Uncompressed 2-bit binary values encoded in an array of integers
 * @param U_len Unumber of integers in U
 * @param used_bits Number of bits in U that are used
 * @param B Resulting bitmap index where the index for value 0 starts at
 *          position 0 and the index for value 1 starts at position N, value 2
 *          at N*2, etc. Where U encodeds N (U_len*16) 2-bit values
 * 
 * @retval size total number of elements in B, where the bitmap index for each
 *         value uses 1/4 of the size
 *
 * Example Usage:
 * @code
 *      uint32_t U1[4] = {
 *              bin_char_to_int("00011011000110110001101100011011"),
 *              bin_char_to_int("00000101010110101111000001010110"),
 *              bin_char_to_int("00000000000000000000000000000000"),
 *              bin_char_to_int("11111111111111111111111111111111")
 *      };
 *      uint32_t *B;
 *      uint32_t B_len = ubin_to_bitmap(U1,4,128,u,&B);
 * @endcode
 */
uint32_t  ubin_to_bitmap(uint32_t *U,
                             uint32_t U_len,
                             uint32_t used_bits,
                             uint32_t **B);

/**
 * @brief Convert an array of uncompressed binary values to a bitmap index of
 *        the values using 16-bit ints.
 *
 * @param U Uncompressed 2-bit binary values encoded in an array of integers
 * @param U_len Unumber of integers in U
 * @param used_bits Number of bits in U that are used
 * @param B Resulting bitmap index where the index for value 0 starts at
 *          position 0 and the index for value 1 starts at position N, value 2
 *          at N*2, etc. Where U encodeds N (U_len*16) 2-bit values
 * 
 * @retval size total number of elements in B, where the bitmap index for each
 *         value uses 1/4 of the size
 *
 * Example Usage:
 * @code
 *      #include <stdint.h>
 *      uint32_t U1[4] = {
 *              bin_char_to_int("00011011000110110001101100011011"),
 *              bin_char_to_int("00000101010110101111000001010110"),
 *              bin_char_to_int("00000000000000000000000000000000"),
 *              bin_char_to_int("11111111111111111111111111111111")
 *      };
 *      uint16_t *B;
 *      uint32_t B_len = ubin_to_bitmap(U1,4,128,u,&B);
 * @endcode
 */

uint32_t ubin_to_bitmap_wah16(uint32_t *U,
                                  uint32_t U_len,
                                  uint32_t num_fields,
                                  uint16_t **W,
                                  uint32_t **wah_sizes);


/**
 * @brief Convert an array of uncompressed two-bit binary values to WAH
 *        encoding of the bitmap index representation of the two-bit binary
 *
 *  In this scheme, there are 4 uniq values in the uncompressed binary
 *  (0,1,2,3), so the WAH encoding will create 4 different sets, one for each
 *  bitmap index
 *
 * @param U         Uncompressed 2-bit binary values encoded in an array of
 *                  integers
 * @param U_len     Unumber of integers in U
 * @param num_fields    The number of fields
 * @param W         The resulting array of WAH words (memory allocted in 
 *                  function, value set in function)
 * @param wah_sizes A list of the 4 sizes 
 *
 * @retval number of ints in W
 *
 * Example Usage:
 * @code
 *    char *plt = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
 *                "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
 *                "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
 *                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
 *                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
 *                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
 *                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
 *                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
 *                "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 "
 *                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 "
 *                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 "
 *                "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3";
 *
 *    uint32_t *ubin;
 *    uint32_t ubin_len = plt_line_to_packed_ints(plt, 192, &ubin);
 *
 *    TEST_ASSERT_EQUAL(12, ubin_len);
 *
 *    uint32_t *wah;
 *    uint32_t *wah_sizes;
 *    uint32_t wah_len = ubin_to_bitmap_wah(ubin,
 *                                              ubin_len,
 *                                              192,
 *                                              &wah,
 *                                              &wah_sizes);
 * @endcode
 */
uint32_t ubin_to_bitmap_wah(uint32_t *U,
                            uint32_t U_len,
                            uint32_t num_fields,
                            uint32_t **W,
                            uint32_t **wah_sizes);

/**
 * @brief convert an uncompressed binary file to a WAH encoded bitmap index
 *
 * A WAH file is: 
 * number of fields (32-bit)
 * number of records (32-bit)
 * Record offsets (number of records * 4 * 32-bit)
 * 1st record 0 bitmap
 * 1st record 1 bitmap
 * 1st record 2 bitmap
 * 1st record 3 bitmap
 * 2nd record 0 bitmap
 * 2nd record 1 bitmap
 * 2nd record 2 bitmap
 * 2nd record 3 bitmap
 *
 * @param ubin_in uncompressed binary file name
 * @param wah_out WAH encoded bitmap index file
 *
 * @retval 0 if everything went right
 * @retval 1 otherwise
 *
 * @code
 *      char *plt_file_name="data/10.1e4.ind.txt";
 *      char *ubin_file_name="data/10.1e4.ind.ubin";
 *      char *wah_file_name="data/10.1e4.ind.wah";
 *
 *      convert_file_by_name_plt_to_ubin(plt_file_name, ubin_file_name);
 *      convert_file_by_name_ubin_to_wahbm(ubin_file_name, wah_file_name);
 * @endcode
 */
uint32_t convert_file_by_name_ubin_to_wahbm(char *ubin_in,
                                            char *wah_out,
                                            char *full_cmd);

/**
 * @brief convert an uncompressed binary file to a 16-bit WAH encoded bitmap
 * index
 *
 * A WAH file is: 
 * number of fields (32-bit)
 * number of records (32-bit)
 * Record offsets (number of records * 4 * 16-bit)
 * 1st record 0 bitmap
 * 1st record 1 bitmap
 * 1st record 2 bitmap
 * 1st record 3 bitmap
 * 2nd record 0 bitmap
 * 2nd record 1 bitmap
 * 2nd record 2 bitmap
 * 2nd record 3 bitmap
 *
 * @param ubin_in uncompressed binary file name
 * @param wah_out 16-bit WAH encoded bitmap index file
 *
 * @retval 0 if everything went right
 * @retval 1 otherwise
 *
 * @code
 *      char *plt_file_name="data/10.1e4.ind.txt";
 *      char *ubin_file_name="data/10.1e4.ind.ubin";
 *      char *wah_file_name="data/10.1e4.ind.wah16";
 *
 *      convert_file_by_name_plt_to_ubin(plt_file_name, ubin_file_name);
 *      convert_file_by_name_ubin_to_wahbm(ubin_file_name, wah_file_name);
 * @endcode
 */
uint32_t convert_file_by_name_ubin_to_wahbm16(char *ubin_in,
                                                  char *wah_out);
/**
 * @brief Get a pointer to the uncompressed binary encoded record
 *
 * @param uf The ubin data structure
 * @param record_id The index of the record of interest
 * @param ubin_record A pointer to the uncompressed binary record
 *
 * @retval number of integers in the record
 *
 * Example Usage:
 * @code
 *     char *ubin_file_name="data/10.1e4.ind.ubin";
 *     struct ubin_file uf = init_ubin_file(ubin_file_name);
 *
 *     uint32_t *ints, num_ints;
 *     uint32_t record_id = 0;
 *     num_ints = get_ubin_record(uf, record_id, &ints);
 * @endcode
 */
uint32_t get_ubin_record(struct ubin_file uf,
                             uint32_t record_id,
                             uint32_t **ubin_record);
/**
 * @brief
 *
 * @param
 *
 * @returnval
 *
 * Example Usage:
 * @code
 * @endcode
 */

uint32_t gt_records_ubin(struct ubin_file uf,
                             uint32_t *record_ids,
                             uint32_t num_r,
                             uint32_t test_value,
                             uint32_t **R);

/**
 * @brief convert an uncompressed binary file using WAH (no bitmaps)
 *
 * A WAH file is: 
 * number of fields (32-bit)
 * number of records (32-bit)
 * Record offsets (number of records * 32-bit)
 * 1st record 
 * 2nd record 
 * ...
 *
 * @param ubin_in uncompressed binary file name
 * @param wah_out WAH encoded ubin file
 *
 * @retval 0 if everything went right
 * @retval 1 otherwise
 *
 * @code
 *      char *plt_file_name="data/10.1e4.ind.txt";
 *      char *ubin_file_name="data/10.1e4.ind.ubin";
 *      char *wah_file_name="data/10.1e4.ind.wah";
 *
 *      convert_file_by_name_plt_to_ubin(plt_file_name, ubin_file_name);
 *      convert_file_by_name_ubin_to_wah(ubin_file_name, wah_file_name);
 * @endcode
 */
uint32_t convert_file_by_name_ubin_to_wah(char *ubin_in, char *wah_out);


/**
 * @brief Print an uncompressed binary file.
 *
 * If num_r > 0, then record_ids should contain the ids of records in the file
 * that will be displayed, otherwise all records will be displayed.
 *
 * @param uf An intitialized uncompressed binary file
 * @param record_ids An array of record ids
 * @param num_r number of records in the array
 * @param format Output format: 0:plain text
 *                              1:packed int
 *
 * @retval number of records printed 
 */
uint32_t print_ubin(struct ubin_file uf,
                        uint32_t *record_ids,
                        uint32_t num_r,
                        uint32_t format);


/**
 * @brief Print uncompressed binary file (wrapper around print_ubin). 
 *
 * If num_r > 0, then record_ids should contain the ids of records in the file
 * that will be displayed, otherwise all records will be displayed.
 *
 * @param ubin_file_name Uncompressed binary file name
 * @param record_ids An array of record ids
 * @param num_r number of records in the array
 * @param format Output format: 0:plain text
 *                              1:packed int
 *
 *
 * @retval number of records printed 
 */
uint32_t print_by_name_ubin(char *ubin_file_name,
                               uint32_t *record_ids,
                               uint32_t num_r,
                               uint32_t format);

/**
 * @brief Return records whose values are >= start_test_value and <
 * end_test_value
 *
 * @param uf The initialized uncompressed binary data structure
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
uint32_t range_records_ubin(struct ubin_file uf,
                            uint32_t *record_ids,
                            uint32_t num_r,
                            uint32_t start_test_value,
                            uint32_t end_test_value,
                            uint32_t **R);

uint32_t count_range_records_ubin(struct ubin_file uf,
                                  uint32_t *record_ids,
                                  uint32_t num_r,
                                  uint32_t start_test_value,
                                  uint32_t end_test_value,
                                  uint32_t **R);

uint32_t convert_file_by_name_ubin_to_plt(char *ubin_in, char *plt_out);

/*
uint32_t convert_hdf5_ind_ubin_to_ind_wah(struct hdf5_file hdf5_f,
                                              char *wah_out);
*/
#endif
