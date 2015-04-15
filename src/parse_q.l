%{
#include "parse_q.h"
%}
%option nounput yylineno

%%

"count"             return p_count;
"pct"               return p_pct;
"maf"               return p_maf;
"=="                return p_equal;
"!="                return p_not_equal;
"<"                 return p_less_than;
">"                 return p_greater_than;
"<="                return p_less_than_equal;
">="                return p_greater_than_equal;
"HET"               return p_het;
"HOM_REF"          return p_homo_ref;
"HOM_ALT"          return p_homo_alt;
"UNKNOWN"           return p_unknown;
"("                 return p_r_paren;
")"                 return p_l_paren;
[0-9]+(\.[0-9]*)?   return p_number;
[ \t\n]             ;


%%

int yywrap(void)
{
    return 1;
}
