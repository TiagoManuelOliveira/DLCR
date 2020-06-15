#include "msg.h"
#include <stdio.h>
#include <stdlib.h>

void ModelsExplanation(void){
  fprintf(stderr,
  "                                                                        \n"		  
  "      -m [NB_C]:[NB_D]:[NB_I]:[NB_H]:[NB_G]/[NB_S]:[NB_E]:[NB_A]        \n"
  "           Template of a target context model.                          \n"
  "           Parameters:                                                  \n"
  "           [NB_C]: (integer [1;20]) order size of the regular context   \n"
  "                   model. Higher values use more RAM but, usually, are  \n"
  "                   related to a better compression score.               \n"
  "           [NB_D]: (integer [1;5000]) denominator to build alpha, which \n"
  "                   is a parameter estimator. Alpha is given by 1/[NB_D].\n"
  "                   Higher values are usually used with higher [NB_C],   \n"
  "                   and related to confiant bets. When [NB_D] is one,    \n"
  "                   the probabilities assume a Laplacian distribution.   \n"
  "           [NB_I]: (integer {0,1,2}) number to define if a sub-program  \n"
  "                   which addresses the specific properties of DNA       \n"
  "                   sequences (Inverted repeats) is used or not. The     \n"
  "                   number 2 turns ON this sub-program without the       \n"
  "                   regular context model (only inverted repeats). The   \n"
  "                   number 1 turns ON the sub-program using at the same  \n"
  "                   time the regular context model. The number 0 does    \n"
  "                   not contemple its use (Inverted repeats OFF). The    \n"
  "                   use of this sub-program increases the necessary time \n"
  "                   to compress but it does not affect the RAM.          \n"
  "           [NB_H]: (integer [1;254]) size of the cache-hash for deeper  \n"
  "                   context models, namely for [NB_C] > 14. When the     \n"
  "                   [NB_C] <= 14 use, for example, 1 as a default. The   \n"
  "                   RAM is highly dependent of this value (higher value  \n"
  "                   stand for higher RAM).                               \n"
  "           [NB_G]: (real [0;1)) real number to define gamma. This value \n"
  "                   represents the decayment forgetting factor of the    \n"
  "                   regular context model in definition.                 \n"
  "           [NB_S]: (integer [0;20]) maximum number of editions allowed  \n"
  "                   to use a substitutional tolerant model with the same \n"
  "                   memory model of the regular context model with       \n"
  "                   order size equal to [NB_C]. The value 0 stands for   \n"
  "                   turning the tolerant context model off. When the     \n"
  "                   model is on, it pauses when the number of editions   \n"
  "                   is higher that [NB_C], while it is turned on when    \n"
  "                   a complete match of size [NB_C] is seen again. This  \n"
  "                   is probabilistic-algorithmic model very usefull to   \n"
  "                   handle the high substitutional nature of genomic     \n"
  "                   sequences. When [NB_S] > 0, the compressor used more \n"
  "                   processing time, but uses the same RAM and, usually, \n"
  "                   achieves a substantial higher compression ratio. The \n"
  "                   impact of this model is usually only noticed for     \n"
  "                   [NB_C] >= 14.                                        \n"
  "           [NB_E]: (integer [1;5000]) denominator to build alpha for    \n"
  "                   substitutional tolerant context model. It is         \n"
  "                   analogous to [NB_D], however to be only used in the  \n"
  "                   probabilistic model for computing the statistics of  \n"
  "                   the substitutional tolerant context model.           \n"
  "           [NB_A]: (real [0;1)) real number to define gamma. This value \n"
  "                   represents the decayment forgetting factor of the    \n"
  "                   substitutional tolerant context model in definition. \n"
  "                   Its definition and use is analogus to [NB_G].        \n"
  "                                                                        \n");
  }

void PrintMenuCompression(void){
  fprintf(stderr,
  "                                                                        \n"
  "      DLCR: Efficient detection of Distant Low Complexity Regions       \n"
  "      ===========================================================       \n"
  "                                                                        \n"
  "AUTHORS                                                                 \n"
  "      Tiago Oliveira (tiagomanuel28@gmail.com) and D. Pratas            \n"
  "                                                                        \n"
  "SYNOPSIS                                                                \n"
  "      ./DLCR [OPTION]... [FILE]                                         \n"
  "                                                                        \n"
  "SAMPLE                                                                  \n"
  "      ./DLCR -v -l 3 -r 512 -t 1 -w 0.025 genome.fa                     \n"
  "                                                                        \n"
  "DESCRIPTION                                                             \n"
  "      Quantification and Localization of                                \n"
  "      Distant Low Complexity Regions in FASTA files.                    \n"
  "                                                                        \n"
  "      -h,  --help                                                       \n"
  "           usage guide (help menu).                                     \n"
  "                                                                        \n"
  "      -V,  --version                                                    \n"
  "           Display program and version information.                     \n"
  "                                                                        \n"
  "      -F,  --force                                                      \n"
  "           force mode. Overwrites old files.                            \n"
  "                                                                        \n"
  "      -v,  --verbose                                                    \n"
  "           verbose mode (more information).                             \n"
  "                                                                        \n"
  "      -t [NUMBER],  --threshold [NUMBER]                                \n"
  "           Threshold to segment regions (real).                         \n"
  "                                                                        \n"
  "      -w [NUMBER],  --weight [NUMBER]                                   \n"
  "           Weight to use in low-pass filter (real).                     \n"
  "                                                                        \n"
  "      -i [NUMBER],  --ignore [NUMBER]                                   \n"
  "           Ignore lengths of segmented regions below this value.        \n"
  "                                                                        \n"
  "      -r [NUMBER],  --region-size [NUMBER]                              \n"
  "           Region size to ignore while updating the models.             \n"
  "           The latest region-size is not updated in the models.         \n"
  "                                                                        \n"
  "      -p,  --show-parameters                                            \n"
  "           show parameters of the models for optimization.              \n"
  "                                                                        \n"
  "      -s,  --show-levels                                                \n"
  "           show pre-computed compression levels (configured parameters).\n"
  "                                                                        \n",
  VERSION, RELEASE);

  fprintf(stderr,
  "      -l [NUMBER],  --level [NUMBER]                                    \n"
  "           Compression level (integer).                                 \n"
  "           Default level: %u.                                           \n"
  "           It defines compressibility in balance with computational     \n"
  "           resources (RAM & time). Use -s for levels perception.        \n",
  DEFAULT_LEVEL);

  fprintf(stderr,
  "                                                                        \n"
  "      [FILE]                                                            \n"
  "           Input sequence filename (to analyze) -- MANDATORY.           \n"
  "           File to analyze (last argument).                             \n"
  "                                                                        \n"
  "COPYRIGHT                                                               \n"
  "      Copyright (C) 2020, IEETA, University of Aveiro.                  \n"
  "      This is a Free software, under GPLv3. You may redistribute        \n"
  "      copies of it under the terms of the GNU - General Public          \n"
  "      License v3 <http://www.gnu.org/licenses/gpl.html>. There          \n"
  "      is NOT ANY WARRANTY, to the extent permitted by law.              \n"
  "                                                                        \n");
  }

void PrintVersion(void){
  fprintf(stderr,
  "                                                                        \n"
  "                          =================                             \n"
  "                          |    DLCR %u.%u   |                           \n"
  "                          =================                             \n"
  "                                                                        \n"
  "             Detection of Distant Low Complexity Regions.               \n"
  "                                                                        \n"
  "               Copyright (C) 2020 University of Aveiro.                 \n"
  "                                                                        \n"
  "                This is a Free software, under GPLv3.                   \n"
  "                                                                        \n"
  "You may redistribute copies of it under the terms of the GNU - General  \n"
  "Public License v3 <http://www.gnu.org/licenses/gpl.html>. There is NOT  \n"
  "ANY WARRANTY, to the extent permitted by law. Developed and Written by  \n"
  "Diogo Pratas and Tiago M. Oliveira.\n\n", VERSION, RELEASE);
  }
