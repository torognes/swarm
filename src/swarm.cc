/*
    SWARM

    Copyright (C) 2012-2017 Torbjorn Rognes and Frederic Mahe

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

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#include "swarm.h"

/* ARGUMENTS AND THEIR DEFAULTS */

#define DEFAULT_GAPOPEN 12
#define DEFAULT_GAPEXTEND 4
#define DEFAULT_MATCHSCORE 5
#define DEFAULT_MISMATCHSCORE (-4)
#define DEFAULT_THREADS 1
#define DEFAULT_RESOLUTION 1
#define DEFAULT_BREAK_SWARMS 0
#define DEFAULT_MOTHUR 0
#define DEFAULT_USEARCH_ABUNDANCE 0
#define DEFAULT_INTERNAL_STRUCTURE 0
#define DEFAULT_LOG 0
#define DEFAULT_NO_OTU_BREAKING 0
#define DEFAULT_FASTIDIOUS 0
#define DEFAULT_BOUNDARY 3

char * outfilename;
char * statsfilename;
char * uclustfilename;
char * progname;
char * databasename;
long gapopen;
long gapextend;
long matchscore;
long mismatchscore;
unsigned long threads;
long resolution;
long break_swarms;
long mothur;
long usearch_abundance;

char * opt_log;
char * opt_internal_structure;
char * opt_seeds;
long opt_no_otu_breaking;
long opt_fastidious;
long opt_boundary;
long opt_bloom_bits;
long opt_ceiling;
long opt_append_abundance;

long penalty_factor;
long penalty_gapextend;
long penalty_gapopen;
long penalty_mismatch;

/* Other variables */

long mmx_present = 0;
long sse_present = 0;
long sse2_present = 0;
long sse3_present = 0;
long ssse3_present = 0;
long sse41_present = 0;
long sse42_present = 0;
long popcnt_present = 0;
long avx_present = 0;
long avx2_present = 0;

unsigned long dbsequencecount = 0;

FILE * outfile;
FILE * statsfile;
FILE * uclustfile;
FILE * logfile = stderr;
FILE * internal_structure_file;
FILE * fp_seeds = 0;

char sym_nt[] = "-acgt                           ";

void cpuid(unsigned int f1,
           unsigned int f2,
           unsigned int & a,
           unsigned int & b,
           unsigned int & c,
           unsigned int & d)
{
  __asm__ __volatile__ ("cpuid"
                        : "=a" (a), "=b" (b), "=c" (c), "=d" (d)
                        : "a" (f1), "c" (f2));
}

void cpu_features_detect()
{
  unsigned int a, b, c, d;

  cpuid(0, 0, a, b, c, d);
  unsigned int maxlevel = a & 0xff;

  if (maxlevel >= 1)
  {
    cpuid(1, 0, a, b, c, d);
    mmx_present    = (d >> 23) & 1;
    sse_present    = (d >> 25) & 1;
    sse2_present   = (d >> 26) & 1;
    sse3_present   = (c >>  0) & 1;
    ssse3_present  = (c >>  9) & 1;
    sse41_present  = (c >> 19) & 1;
    sse42_present  = (c >> 20) & 1;
    popcnt_present = (c >> 23) & 1;
    avx_present    = (c >> 28) & 1;

    if (maxlevel >= 7)
    {
      cpuid(7, 0, a, b, c, d);
      avx2_present   = (b >>  5) & 1;
    }
  }
}

void cpu_features_show()
{
  fprintf(logfile, "CPU features:     ");
  if (mmx_present)
    fprintf(logfile, " mmx");
  if (sse_present)
    fprintf(logfile, " sse");
  if (sse2_present)
    fprintf(logfile, " sse2");
  if (sse3_present)
    fprintf(logfile, " sse3");
  if (ssse3_present)
    fprintf(logfile, " ssse3");
  if (sse41_present)
    fprintf(logfile, " sse4.1");
  if (sse42_present)
    fprintf(logfile, " sse4.2");
  if (popcnt_present)
    fprintf(logfile, " popcnt");
  if (avx_present)
    fprintf(logfile, " avx");
  if (avx2_present)
    fprintf(logfile, " avx2");
  fprintf(logfile, "\n");
}

long args_long(char * str, const char * option)
{
  char * endptr;
  long temp = strtol(str, & endptr, 10);
  if (*endptr)
    fatal("Invalid numeric argument for option %s", option);
  return temp;
}

void args_show()
{
  cpu_features_show();
  fprintf(logfile, "Database file:     %s\n", databasename ? databasename : "(stdin)");
  fprintf(logfile, "Output file:       %s\n", outfilename ? outfilename : "(stdout)");
  if (statsfilename)
    fprintf(logfile, "Statistics file:   %s\n", statsfilename);
  if (uclustfilename)
    fprintf(logfile, "Uclust file:       %s\n", uclustfilename);
  fprintf(logfile, "Resolution (d):    %ld\n", resolution);
  fprintf(logfile, "Threads:           %lu\n", threads);

  if (resolution > 1)
    {
  fprintf(logfile, "Scores:            match: %ld, mismatch: %ld\n", matchscore, mismatchscore);
  fprintf(logfile, "Gap penalties:     opening: %ld, extension: %ld\n", gapopen, gapextend);
  fprintf(logfile, "Converted costs:   mismatch: %ld, gap opening: %ld, gap extension: %ld\n", penalty_mismatch, penalty_gapopen, penalty_gapextend);
    }
  fprintf(logfile, "Break OTUs:        %s\n", opt_no_otu_breaking ? "No" : "Yes");
  if (opt_fastidious)
    fprintf(logfile, "Fastidious:        Yes, with boundary %ld\n", opt_boundary);
  else
    fprintf(logfile, "Fastidious:        No\n");
}

void args_usage()
{
  /*               0         1         2         3         4         5         6         7          */
  /*               01234567890123456789012345678901234567890123456789012345678901234567890123456789 */

  fprintf(stderr, "Usage: swarm [OPTIONS] [filename]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "General options:\n");
  fprintf(stderr, " -h, --help                          display this help and exit\n");
  fprintf(stderr, " -t, --threads INTEGER               number of threads to use (1)\n");
  fprintf(stderr, " -v, --version                       display version information and exit\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Clustering options:\n");
  fprintf(stderr, " -b, --boundary INTEGER              min mass of large OTU for fastidious (3)\n");
  fprintf(stderr, " -c, --ceiling INTEGER               max memory in MB used for fastidious\n");
  fprintf(stderr, " -d, --differences INTEGER           resolution (1)\n");
  fprintf(stderr, " -f, --fastidious                    link nearby low-abundance swarms\n");
  fprintf(stderr, " -n, --no-otu-breaking               never break OTUs\n");
  fprintf(stderr, " -y, --bloom-bits INTEGER            bits used per Bloom filter entry (16)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Input/output options:\n");
  fprintf(stderr, " -a, --append-abundance INTEGER      value to use when abundance is missing\n");
  fprintf(stderr, " -i, --internal-structure FILENAME   write internal swarm structure to file\n");
  fprintf(stderr, " -l, --log FILENAME                  log to file, not to stderr\n");
  fprintf(stderr, " -o, --output-file FILENAME          output result filename (stdout)\n");
  fprintf(stderr, " -r, --mothur                        output in mothur list file format\n");
  fprintf(stderr, " -s, --statistics-file FILENAME      dump OTU statistics to file\n");
  fprintf(stderr, " -u, --uclust-file FILENAME          output in UCLUST-like format to file\n");
  fprintf(stderr, " -w, --seeds FILENAME                write seed seqs with abundances to FASTA\n");
  fprintf(stderr, " -z, --usearch-abundance             abundance annotation in usearch style\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Pairwise alignment advanced options:\n");
  fprintf(stderr, " -m, --match-reward INTEGER          reward for nucleotide match (5)\n");
  fprintf(stderr, " -p, --mismatch-penalty INTEGER      penalty for nucleotide mismatch (4)\n");
  fprintf(stderr, " -g, --gap-opening-penalty INTEGER   gap open penalty (12)\n");
  fprintf(stderr, " -e, --gap-extension-penalty INTEGER gap extension penalty (4)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "See 'man swarm' for more details.\n");
}

void show_header()
{
  char title[] = "Swarm " SWARM_VERSION;
  char ref[] = "Copyright (C) 2012-2017 Torbjorn Rognes and Frederic Mahe";
  char url[] = "https://github.com/torognes/swarm";
  fprintf(logfile, "%s [%s %s]\n%s\n%s\n\n",
          title, __DATE__, __TIME__, ref, url);
  fprintf(logfile, "Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M (2014)\n");
  fprintf(logfile, "Swarm: robust and fast clustering method for amplicon-based studies\n");
  fprintf(logfile, "PeerJ 2:e593 https://doi.org/10.7717/peerj.593\n");
  fprintf(logfile, "\n");
  fprintf(logfile, "Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M (2015)\n");
  fprintf(logfile, "Swarm v2: highly-scalable and high-resolution amplicon clustering\n");
  fprintf(logfile, "PeerJ 3:e1420 https://doi.org/10.7717/peerj.1420\n");
  fprintf(logfile, "\n");
}

void args_init(int argc, char **argv)
{
  /* Set defaults */

  progname = argv[0];

  databasename = NULL;
  outfilename = NULL;
  statsfilename = NULL;
  resolution = DEFAULT_RESOLUTION;
  threads = DEFAULT_THREADS;
  matchscore = DEFAULT_MATCHSCORE;
  mismatchscore = DEFAULT_MISMATCHSCORE;
  gapopen = DEFAULT_GAPOPEN;
  gapextend = DEFAULT_GAPEXTEND;
  mothur = DEFAULT_MOTHUR;
  usearch_abundance = DEFAULT_USEARCH_ABUNDANCE;
  opt_log = DEFAULT_LOG;
  opt_internal_structure = DEFAULT_INTERNAL_STRUCTURE;
  opt_no_otu_breaking = DEFAULT_NO_OTU_BREAKING;
  opt_fastidious = DEFAULT_FASTIDIOUS;
  opt_boundary = DEFAULT_BOUNDARY;
  opt_bloom_bits = 16;
  opt_seeds = 0;
  opt_ceiling = 0;
  opt_append_abundance = 0;

  opterr = 1;

  char short_options[] = "d:ho:t:vm:p:g:e:s:u:rzi:l:nfb:w:y:c:a:";

  /* unused short option letters: jkqx */

  static struct option long_options[] =
  {
    {"differences",           required_argument, NULL, 'd' },
    {"help",                  no_argument,       NULL, 'h' },
    {"output-file",           required_argument, NULL, 'o' },
    {"threads",               required_argument, NULL, 't' },
    {"version",               no_argument,       NULL, 'v' },
    {"match-reward",          required_argument, NULL, 'm' },
    {"mismatch-penalty",      required_argument, NULL, 'p' },
    {"gap-opening-penalty",   required_argument, NULL, 'g' },
    {"gap-extension-penalty", required_argument, NULL, 'e' },
    {"statistics-file",       required_argument, NULL, 's' },
    {"uclust-file",           required_argument, NULL, 'u' },
    {"mothur",                no_argument,       NULL, 'r' },
    {"usearch-abundance",     no_argument,       NULL, 'z' },
    {"internal-structure",    required_argument, NULL, 'i' },
    {"log",                   required_argument, NULL, 'l' },
    {"no-otu-breaking",       no_argument,       NULL, 'n' },
    {"fastidious",            no_argument,       NULL, 'f' },
    {"boundary",              required_argument, NULL, 'b' },
    {"seeds",                 required_argument, NULL, 'w' },
    {"bloom-bits",            required_argument, NULL, 'y' },
    {"ceiling",               required_argument, NULL, 'c' },
    {"append-abundance",      required_argument, NULL, 'a' },
    { 0, 0, 0, 0 }
  };
  
  int used_options[26] = { 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0 };

  int option_index = 0;
  int c;
  
  while ((c = getopt_long(argc, argv, short_options, long_options, &option_index)) != -1)
  {
    if ((c >= 'a') && (c <= 'z'))
      {
        int optindex = c - 'a';
        if (used_options[optindex] == 1)
          {
            int longoptindex = 0;
            while (long_options[longoptindex].name)
              {
                if (long_options[longoptindex].val == c)
                  break;
                longoptindex++;
              }

            fprintf(stderr,
                    "WARNING: Option -%c or --%s specified more than once.\n",
                    c,
                    long_options[longoptindex].name);
          }
        used_options[optindex]++;
      }

    switch(c)
    {
    case 'd':
      /* differences (resolution) */
      resolution = args_long(optarg, "-d or --differences");
      break;
          
    case 'h':
      /* help */
      show_header();
      args_usage();
      exit(0);
      break;

    case 'o':
      /* output-file */
      outfilename = optarg;
      break;
          
    case 't':
      /* threads */
      threads = args_long(optarg, "-t or --threads");
      break;
          
    case 'v':
      /* version */
      show_header();
      exit(0);
      break;

    case 'm':
      /* match-reward */
      matchscore = args_long(optarg, "-m or --match-reward");
      break;
          
    case 'p':
      /* mismatch-penalty */
      mismatchscore = - args_long(optarg, "-p or --mismatch-penalty");
      break;
          
    case 'g':
      /* gap-opening-penalty */
      gapopen = args_long(optarg, "-g or --gap-opening-penalty");
      break;
          
    case 'e':
      /* gap extension penalty */
      gapextend = args_long(optarg, "-e or --gap-extension-penalty");
      break;
          
    case 's':
      /* statistics-file */
      statsfilename = optarg;
      break;
          
    case 'u':
      /* uclust-file */
      uclustfilename = optarg;
      break;
          
    case 'r':
      /* mothur */
      mothur = 1;
      break;
          
    case 'z':
      /* usearch-abundance */
      usearch_abundance = 1;
      break;
          
    case 'i':
      /* internal-structure */
      opt_internal_structure = optarg;
      break;
          
    case 'l':
      /* log */
      opt_log = optarg;
      break;
          
    case 'n':
      /* no-otu-breaking */
      opt_no_otu_breaking = 1;
      break;
          
    case 'f':
      /* fastidious */
      opt_fastidious = 1;
      break;
          
    case 'b':
      /* boundary */
      opt_boundary = args_long(optarg, "-b or --boundary");
      break;
          
    case 'w':
      /* seeds */
      opt_seeds = optarg;
      break;
      
    case 'y':
      /* bloom-bits */
      opt_bloom_bits = args_long(optarg, "-y or --bloom-bits");
      break;
      
    case 'c':
      /* ceiling */
      opt_ceiling = args_long(optarg, "-c or --ceiling");
      break;
      
    case 'a':
      /* append-abundance */
      opt_append_abundance = args_long(optarg, "-a or --append-abundance");
      break;
      
    default:
      show_header();
      args_usage();
      exit(1);
      break;
    }
  }
  
  if (optind < argc)
    databasename = argv[optind];
  
  if ((resolution < 0) || (resolution > 255))
    fatal("Error: number of differences specified with -d must be in the range 0 to 255.");

  if ((threads < 1) || (threads > MAX_THREADS))
    fatal("Illegal number of threads specified");
  
  if ((gapopen < 0) || (gapextend < 0) || ((gapopen + gapextend) < 1))
    fatal("Illegal gap penalties specified.");

  if (matchscore < 1)
    fatal("Illegal match reward specified.");

  if (mismatchscore > -1)
    fatal("Illegal mismatch penalty specified.");

  if ((opt_bloom_bits < 2) || (opt_bloom_bits > 64))
    fatal("Illegal number of Bloom filter bits specified (must be 2..64).");

  if (opt_ceiling < 0)
    fatal("Illegal memory ceiling");

  if (opt_append_abundance < 0)
    fatal("Illegal abundance value specified");

  if ((opt_ceiling > 0) && (opt_fastidious == 0))
    fprintf(stderr, "WARNING: Options -c and --ceiling ignored without -f or --fastidious.\n");

  if (outfilename)
    {
      outfile = fopen(outfilename, "w");
      if (! outfile)
        fatal("Unable to open output file for writing.");
    }
  else
    outfile = stdout;
  
  if (opt_seeds)
    {
      fp_seeds = fopen(opt_seeds, "w");
      if (! fp_seeds)
        fatal("Unable to open seeds file for writing.");
    }
  else
    fp_seeds = 0;
  
  if (statsfilename)
    {
      statsfile = fopen(statsfilename, "w");
      if (! statsfile)
        fatal("Unable to open statistics file for writing.");
    }
  else
    statsfile = 0;
  
  if (uclustfilename)
    {
      uclustfile = fopen(uclustfilename, "w");
      if (! uclustfile)
        fatal("Unable to open uclust file for writing.");
    }
  else
    uclustfile = 0;

  if (opt_log)
    {
      logfile = fopen(opt_log, "w");
      if (! logfile)
        fatal("Unable to open log file for writing.");
    }
  else
    logfile = stderr;

  if (opt_internal_structure)
    {
      internal_structure_file = fopen(opt_internal_structure, "w");
      if (! internal_structure_file)
        fatal("Unable to open internal structure file for writing.");
    }
  else
    internal_structure_file = stderr;

  if (opt_fastidious && (resolution != 1))
    fatal("The fastidious option only works when the resolution (d) is 1.\n");
}

int main(int argc, char** argv)
{
  cpu_features_detect();

  if (!sse2_present)
    fatal("This program requires a processor with SSE2 instructions.\n");
  
  args_init(argc, argv);

  penalty_mismatch = 2 * matchscore - 2 * mismatchscore;
  penalty_gapopen = 2 * gapopen;
  penalty_gapextend = matchscore + 2 * gapextend;

  penalty_factor = gcd(gcd(penalty_mismatch, penalty_gapopen), penalty_gapextend);
  
  penalty_mismatch /= penalty_factor;
  penalty_gapopen /= penalty_factor;
  penalty_gapextend /= penalty_factor;

  show_header();
  
  args_show();

  fprintf(logfile, "\n");

  db_read(databasename);
  
  fprintf(logfile, "Database info:     %lu nt", db_getnucleotidecount());
  fprintf(logfile, " in %lu sequences,", db_getsequencecount());
  fprintf(logfile, " longest %lu nt\n", db_getlongestsequence());

  dbsequencecount = db_getsequencecount();

  score_matrix_init();

  search_begin();
  
  switch (resolution)
    {
    case 0:
      dereplicate();
      break;
      
    case 1:
      algo_d1_run();
      break;

    default:
      algo_run();
      break;
    }

  search_end();

  score_matrix_free();

  db_free();

  if (opt_seeds)
    fclose(fp_seeds);

  if (uclustfile)
    fclose(uclustfile);

  if (statsfile)
    fclose(statsfile);

  if (outfilename)
    fclose(outfile);
}
