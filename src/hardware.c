/*
    Copyright (C) 2017 Tomas Flouri

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
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "pll.h"

#ifndef __PPC__
#define cpuid(f1, f2, a, b, c, d)                                \
  __asm__ __volatile__ ("cpuid"                                  \
                        : "=a" (a), "=b" (b), "=c" (c), "=d" (d) \
                        : "a" (f1), "c" (f2));
#endif

static void cpu_features_detect()
{
  unsigned int a,b,c,d;

  memset(pll_hardware,0,sizeof(pll_hardware_t));

#ifdef __PPC__
  pll_hardware->altivec_present = 1;
#else

  cpuid(0,0,a,b,c,d);
  unsigned int maxlevel = a & 0xff;

  if (maxlevel >= 1)
  {
    cpuid(1,0,a,b,c,d);
    pll_hardware->mmx_present    = (d >> 23) & 1;
    pll_hardware->sse_present    = (d >> 25) & 1;
    pll_hardware->sse2_present   = (d >> 26) & 1;
    pll_hardware->sse3_present   = (c >>  0) & 1;
    pll_hardware->ssse3_present  = (c >>  9) & 1;
    pll_hardware->sse41_present  = (c >> 19) & 1;
    pll_hardware->sse42_present  = (c >> 20) & 1;
    pll_hardware->popcnt_present = (c >> 23) & 1;
    pll_hardware->avx_present    = (c >> 28) & 1;

    if (maxlevel >= 7)
    {
      cpuid(7,0,a,b,c,d);
      pll_hardware->avx2_present = (b >> 5) & 1;
    }
  }
#endif
}

static void cpu_features_show()
{
  if (!pll_hardware)
  {
    /* TODO: Add proper error control after we figure out
       cross-platform compatibility */
    return;
  }
    
  fprintf(stderr, "CPU features:");
  if (pll_hardware->altivec_present)
    fprintf(stderr, " altivec");
  if (pll_hardware->mmx_present)
    fprintf(stderr, " mmx");
  if (pll_hardware->sse_present)
    fprintf(stderr, " sse");
  if (pll_hardware->sse2_present)
    fprintf(stderr, " sse2");
  if (pll_hardware->sse3_present)
    fprintf(stderr, " sse3");
  if (pll_hardware->ssse3_present)
    fprintf(stderr, " ssse3");
  if (pll_hardware->sse41_present)
    fprintf(stderr, " sse4.1");
  if (pll_hardware->sse42_present)
    fprintf(stderr, " sse4.2");
  if (pll_hardware->popcnt_present)
    fprintf(stderr, " popcnt");
  if (pll_hardware->avx_present)
    fprintf(stderr, " avx");
  if (pll_hardware->avx2_present)
    fprintf(stderr, " avx2");
  fprintf(stderr, "\n");
}

PLL_EXPORT int pll_hardware_probe()
{
  /* probe cpu features */
  if (!pll_hardware)
  {
    if (!(pll_hardware = (pll_hardware_t *)calloc(1,sizeof(pll_hardware_t))))
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
      return PLL_FAILURE;
    }
  }
  cpu_features_detect();

  return PLL_SUCCESS;
}

PLL_EXPORT void pll_hardware_dump()
{
  cpu_features_show(); 
}
