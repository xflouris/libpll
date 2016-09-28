/*
    Copyright (C) 2015 Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/
%{
#include "pll.h"
#include "parse_phylip.h"

static size_t string_length = 0;

static char * xstrndup(const char * s, size_t len)
{
  char * p = (char *)malloc(len+1);
  if (!p)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }
  strncpy(p,s,len);
  p[len] = 0;
  return p;
}

static char * append(size_t * dstlen, const char * src, size_t srclen)
{
  char * mem = (char *)realloc(pll_phylip_lval.s, (*dstlen + srclen + 1));
  if (!mem)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }
  pll_phylip_lval.s =  mem;
  strncpy(&pll_phylip_lval.s[*dstlen], src, srclen);
  pll_phylip_lval.s[*dstlen + srclen] = 0;
  *dstlen += srclen;
  return pll_phylip_lval.s;
}

%}
%option noyywrap
%option prefix="pll_phylip_"
%option nounput
%option noinput
%s apos
%s quot

%%
<quot>{
\\\"                                         append(&string_length, "\"", 1);
\'                                           append(&string_length, "\'", 1);
\"                                           {BEGIN(INITIAL);return STRING;}
}

<apos>{
\\\'                                         append(&string_length, "\'", 1);
\"                                           append(&string_length, "\"", 1);
\'                                           {BEGIN(INITIAL);return STRING;}
}

<apos,quot>{
\\n                                          append(&string_length, "\n", 1);
\\t                                          append(&string_length, "\t", 1);
\\\\                                         append(&string_length, "\\", 1);
([^\"\'\\]|\n)+                              append(&string_length, pll_phylip_text, pll_phylip_leng);
}

\"                                           {string_length =0; pll_phylip_lval.s = NULL; ;BEGIN(quot);}
\'                                           {string_length =0; pll_phylip_lval.s = NULL; ;BEGIN(apos);}
([^\t\n\r ])+                {pll_phylip_lval.s = xstrndup(pll_phylip_text, pll_phylip_leng); return STRING;}
[ \t\n\r]*\#[^\n]*                           ;
[ \t\n\r]                                    ;
.                                            {snprintf(pll_errmsg, 200, "Syntax error (%c)\n", pll_phylip_text[0]);
                                              pll_errno = PLL_ERROR_PHYLIP_SYNTAX;
                                              return PLL_FAILURE;}
%%
