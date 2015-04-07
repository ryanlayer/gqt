#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parse_q.h"

extern int yylex();
extern int yylineo;
extern char *yytext;

typedef struct yy_buffer_state * YY_BUFFER_STATE;
extern int yyparse();
extern YY_BUFFER_STATE yy_scan_string(char * str);
extern void yy_delete_buffer(YY_BUFFER_STATE buffer);

char *names[] = {NULL,
                 "count",
                 "pct",
                 "==",
                 "!=",
                 "<",
                 ">",
                 "<=",
                 ">=",
                 "HET",
                 "HOMO_REF",
                 "HOMO_ALT",
                 "UNKNOWN",
                 "(",
                 ")",
                 "p_number",
                 "maf"};

//{{{ int set_gt(int ntoken, struct gqt_query *q, int *state)
int set_gt(int ntoken, struct gqt_query *q, int *state)
{
    if ((*state & OP_SET_END) ||
        (*state & GT_SET_END) ) {
        fprintf(stderr, "SYNTAX ERROR: invalid placement of '%s' ",
                names[ntoken]);
        return 1;
    }

    switch (ntoken) {
        case p_homo_ref:
            q->genotype_condition[0] = 1;
            return 0;
            break;
        case p_het:
            q->genotype_condition[1] = 1;
            return 0;
            break;
        case p_homo_alt:
            q->genotype_condition[2] = 1;
            return 0;
            break;
        case p_unknown:
            q->genotype_condition[3] = 1;
            return 0;
            break;
        default:
            printf("Syntax error in line\n");
            return 1;
    }
}
//}}}

//{{{ int set_op(int ntoken, struct gqt_query *q, int *state)
int set_op(int ntoken, struct gqt_query *q, int *state)
{
    if (*state != START_OF_EXP) {
        fprintf(stderr, "SYNTAX ERROR: invalid placement of '%s' ",
                names[ntoken]);
        return 1;
    }

    *state |= OP_SET_START;

    q->variant_op = ntoken;

    return 0;
}
//}}}

//{{{ int set_op_cond(int ntoken, struct gqt_query *q, int *state)
int set_op_cond(int ntoken, struct gqt_query *q, int *state)
{
    if ((*state & OP_SET_END) != OP_SET_END) {
        fprintf(stderr, "SYNTAX ERROR: "
                        "Opperation (count,pct,maf) expected prior to '%s' ",
                        names[ntoken]);
        return 1;
    }

    *state |= COND_SET;

    q->op_condition = ntoken;

    return 0;
}
//}}}

//{{{int set_cond_value(char *yytext, struct gqt_query *q, int *state)
int set_cond_value(char *yytext, struct gqt_query *q, int *state)
{
    if ((*state & COND_SET) != COND_SET) {
        fprintf(stderr, "SYNTAX ERROR:"
                        "Opperation (count,pct,maf) and condition "
                        "(==,!=,<, etc.) expected prior to '%s' ",
                        yytext);
        return 1;
    }

    *state |= COND_VALUE_SET;

    q->condition_value = atof(yytext);

    return 0;
}
//}}}

//{{{ int parse_q(char *q_text, struct gqt_query *q_props)
int parse_q(char *q_text, struct gqt_query *q_props)
{
    int ntoken, vtoken;

    int state = START_OF_EXP;

    q_props->variant_op = -1;
    q_props->op_condition = -1;
    q_props->condition_value = -1;

    memset(q_props->genotype_condition,0,4*sizeof(int));

    char *in = (char *) malloc(100*sizeof(char));
    strcpy(in, q_text);
    YY_BUFFER_STATE buffer = yy_scan_string(in);

    ntoken = yylex();

    while (ntoken) {
        switch (ntoken) {
            case p_maf:
            case p_count:
            case p_pct:
                if (set_op(ntoken, q_props, &state))
                    return 1;
                int last_ntoken = ntoken;
                ntoken = yylex();
                if (ntoken != p_r_paren) {
                    fprintf(stderr, "SYNTAX ERROR: '(' expected after '%s' ",
                            names[last_ntoken]);
                    return 1;
                }
                state |= OP_SET_START;
                break;
            case p_het:
            case p_homo_ref:
            case p_homo_alt:
            case p_unknown:
                if (set_gt(ntoken, q_props, &state))
                    return 1;
                state |= GT_SET_START;
                break;
            case p_r_paren:
                fprintf(stderr, "SYNTAX ERROR: "
                            "Opperation (count,pct,maf) expected "
                            "prior to '%s' ",
                            names[ntoken]);
                    return 1;
            case p_l_paren:



                if ( ((q_props->variant_op == p_count) ||
                      (q_props->variant_op == p_pct))  &&
                     (state != (OP_SET_START | GT_SET_START) ) ) {
                    fprintf(stderr, "SYNTAX ERROR: "
                            "Opperation (count,pct) and "
                            "genotype (HOMO_REF,HET,HOMO_ALT,UNKNONW) expected "
                            "prior to '%s' ",
                            names[ntoken]);
                    return 1;
                } else if ( (q_props->variant_op == p_maf) &&
                            (state != OP_SET_START ) ) {
                    fprintf(stderr, "SYNTAX ERROR: "
                            "Opperation (maf) does not expect "
                            "genotype (HOMO_REF,HET,HOMO_ALT,UNKNONW) "
                            "prior to '%s' ",
                            names[ntoken]);
                    return 1;
                }

                state |= OP_SET_END;
                state |= GT_SET_END;
                break;
            case p_equal:
            case p_not_equal:
            case p_less_than:
            case p_greater_than:
            case p_less_than_equal:
            case p_greater_than_equal:
                if (set_op_cond(ntoken, q_props, &state))
                    return 1;
                break;
            case p_number:
                if (set_cond_value(yytext, q_props, &state))
                    return 1;
                break;
            default:
                printf("Syntax error in line\n");
        }

        ntoken = yylex();
    }
    yy_delete_buffer(buffer);

    if (((state & OP_SET_START) == OP_SET_START) &&
        ((state & OP_SET_END) != OP_SET_END)) {
        fprintf(stderr, "SYNTAX ERROR: Missing ')' ");
        return 1;
    }

    if (((state & COND_SET) == COND_SET) &&
        ((state & COND_VALUE_SET) != COND_VALUE_SET)) {
        fprintf(stderr, "SYNTAX ERROR: Missing condition value "
                "(after ==, <, etc.) ");
        return 1;
    }


    /*
    printf("%d\n", q0.variant_op);
    int i;
    for (i = 0; i < 4; ++i)
        printf("%d", q0.genotype_condition[i]);
    printf("\n");
    printf("%d\n", q0.op_condition);
    printf("%f\n", q0.condition_value);
    */

    //free(in);
    return 0;
}
//}}}
