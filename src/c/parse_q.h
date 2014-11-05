#ifndef __PARSE_Q_H__
#define __PARSE_Q_H__

#define     p_count 1
#define     p_pct 2
#define     p_equal 3
#define     p_not_equal 4
#define     p_less_than 5
#define     p_greater_than 6
#define     p_less_than_equal 7
#define     p_greater_than_equal 8
#define     p_het 9
#define     p_homo_ref 10
#define     p_homo_alt 11
#define     p_unknown 12
#define     p_r_paren 13
#define     p_l_paren 14
#define     p_number 15


/* state tracks the expression 
 * From least to most significant
 * 000000 <- start of expression
 * 000001 <- op set start
 * 000010 <- op set end
 * 000100 <- gt set start
 * 001000 <- gt set end
 * 010000 <- cond set
 * 100000 <- cond value set
 */
#define START_OF_EXP    0
#define OP_SET_START    1<<1
#define OP_SET_END      1<<2
#define GT_SET_START    1<<3
#define GT_SET_END      1<<4
#define COND_SET        1<<5
#define COND_VALUE_SET  1<<6

struct gqt_query
{
    int variant_op;
    int genotype_condition[4];
    int op_condition;
    float condition_value;
};

int set_gt(int ntoken, struct gqt_query *q, int *state);

int set_op(int ntoken, struct gqt_query *q, int *state);

int set_op_cond(int ntoken, struct gqt_query *q, int *state);

int parse_q(char *q_text, struct gqt_query *q_props);
#endif
