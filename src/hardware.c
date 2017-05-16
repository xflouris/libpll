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

static void cpu_features_detect()
{
  memset(&pll_hardware,0,sizeof(pll_hardware_t));

  pll_hardware.init            = 1;
#if defined(__PPC__)
  pll_hardware.altivec_present = __builtin_cpu_supports("altivec");
#elif defined(__x86_64__) || defined(__i386__)
  pll_hardware.mmx_present     = __builtin_cpu_supports("mmx");
  pll_hardware.sse_present     = __builtin_cpu_supports("sse");
  pll_hardware.sse2_present    = __builtin_cpu_supports("sse2");
  pll_hardware.sse3_present    = __builtin_cpu_supports("sse3");
  pll_hardware.ssse3_present   = __builtin_cpu_supports("ssse3");
  pll_hardware.sse41_present   = __builtin_cpu_supports("sse4.1");
  pll_hardware.sse42_present   = __builtin_cpu_supports("sse4.2");
  pll_hardware.popcnt_present  = __builtin_cpu_supports("popcnt");
  pll_hardware.avx_present     = __builtin_cpu_supports("avx");
  pll_hardware.avx2_present    = __builtin_cpu_supports("avx2");
#endif
}

static void cpu_features_show()
{
  fprintf(stderr, "CPU features:");
  if (pll_hardware.altivec_present)
    fprintf(stderr, " altivec");
  if (pll_hardware.mmx_present)
    fprintf(stderr, " mmx");
  if (pll_hardware.sse_present)
    fprintf(stderr, " sse");
  if (pll_hardware.sse2_present)
    fprintf(stderr, " sse2");
  if (pll_hardware.sse3_present)
    fprintf(stderr, " sse3");
  if (pll_hardware.ssse3_present)
    fprintf(stderr, " ssse3");
  if (pll_hardware.sse41_present)
    fprintf(stderr, " sse4.1");
  if (pll_hardware.sse42_present)
    fprintf(stderr, " sse4.2");
  if (pll_hardware.popcnt_present)
    fprintf(stderr, " popcnt");
  if (pll_hardware.avx_present)
    fprintf(stderr, " avx");
  if (pll_hardware.avx2_present)
    fprintf(stderr, " avx2");
  fprintf(stderr, "\n");
}

PLL_EXPORT int pll_hardware_probe()
{
  /* probe cpu features */
  cpu_features_detect();

  return PLL_SUCCESS;
}

PLL_EXPORT void pll_hardware_dump()
{
  if (!pll_hardware.init)
    pll_hardware_probe();

  cpu_features_show(); 
}

PLL_EXPORT void pll_hardware_ignore()
{
  pll_hardware.init            = 1;
  pll_hardware.altivec_present = 1;
  pll_hardware.mmx_present     = 1;
  pll_hardware.sse_present     = 1;
  pll_hardware.sse2_present    = 1;
  pll_hardware.sse3_present    = 1;
  pll_hardware.ssse3_present   = 1;
  pll_hardware.sse41_present   = 1;
  pll_hardware.sse42_present   = 1;
  pll_hardware.popcnt_present  = 1;
  pll_hardware.avx_present     = 1;
  pll_hardware.avx2_present    = 1;
}
