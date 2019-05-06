/*
    SWARM

    Copyright (C) 2012-2019 Torbjorn Rognes and Frederic Mahe

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

/* OPTIONS */

char * progname;
char * input_filename;

char * opt_internal_structure;
char * opt_log;
char * opt_output_file;
char * opt_seeds;
char * opt_statistics_file;
char * opt_uclust_file;

long opt_append_abundance;
long opt_bloom_bits;
long opt_boundary;
long opt_ceiling;
long opt_differences;
long opt_fastidious;
long opt_gap_extension_penalty;
long opt_gap_opening_penalty;
long opt_help;
long opt_match_reward;
long opt_mismatch_penalty;
long opt_mothur;
long opt_no_otu_breaking;
long opt_threads;
long opt_usearch_abundance;
long opt_version;

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

unsigned long duplicates_found = 0;

FILE * outfile;
FILE * statsfile;
FILE * uclustfile;
FILE * logfile = stderr;
FILE * internal_structure_file;
FILE * fp_seeds = 0;

char sym_nt[] = "-ACGT                           ";

char * DASH_FILENAME = (char*) "-";
char * STDIN_NAME = (char*) "/dev/stdin";
char * STDOUT_NAME = (char*) "/dev/stdout";

#ifdef __x86_64__

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

#endif

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
#ifdef __x86_64__
  cpu_features_show();
#endif

  fprintf(logfile, "Database file:     %s\n", input_filename);
  fprintf(logfile, "Output file:       %s\n", opt_output_file);
  if (opt_statistics_file)
    fprintf(logfile, "Statistics file:   %s\n", opt_statistics_file);
  if (opt_uclust_file)
    fprintf(logfile, "Uclust file:       %s\n", opt_uclust_file);
  fprintf(logfile, "Resolution (d):    %ld\n", opt_differences);
  fprintf(logfile, "Threads:           %ld\n", opt_threads);

  if (opt_differences > 1)
    {
      fprintf(logfile,
              "Scores:            match: %ld, mismatch: %ld\n",
              opt_match_reward, opt_mismatch_penalty);
      fprintf(logfile,
              "Gap penalties:     opening: %ld, extension: %ld\n",
              opt_gap_opening_penalty, opt_gap_extension_penalty);
      fprintf(logfile,
              "Converted costs:   mismatch: %ld, gap opening: %ld, "
              "gap extension: %ld\n",
              penalty_mismatch, penalty_gapopen, penalty_gapextend);
    }
  fprintf(logfile, "Break OTUs:        %s\n",
          opt_no_otu_breaking ? "No" : "Yes");
  if (opt_fastidious)
    fprintf(logfile, "Fastidious:        Yes, with boundary %ld\n",
            opt_boundary);
  else
    fprintf(logfile, "Fastidious:        No\n");
}

void args_usage()
{
  /*               0         1         2         3         4         5         6         7          */
  /*               01234567890123456789012345678901234567890123456789012345678901234567890123456789 */

  fprintf(stderr, "Usage: swarm [OPTIONS] [FASTAFILE]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "General options:\n");
  fprintf(stderr, " -h, --help                          display this help and exit\n");
  fprintf(stderr, " -t, --threads INTEGER               number of threads to use (1)\n");
  fprintf(stderr, " -v, --version                       display version information and exit\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Clustering options:\n");
  fprintf(stderr, " -d, --differences INTEGER           resolution (1)\n");
  fprintf(stderr, " -n, --no-otu-breaking               never break OTUs (not recommended!)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Fastidious options (only when d = 1):\n");
  fprintf(stderr, " -b, --boundary INTEGER              min mass of large OTU for fastidious (3)\n");
  fprintf(stderr, " -c, --ceiling INTEGER               max memory in MB used for fastidious\n");
  fprintf(stderr, " -f, --fastidious                    link nearby low-abundance swarms\n");
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
  fprintf(stderr, "Pairwise alignment advanced options (only when d > 1):\n");
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
  char ref[] = "Copyright (C) 2012-2019 Torbjorn Rognes and Frederic Mahe";
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
  input_filename = DASH_FILENAME;

  opt_append_abundance = 0;
  opt_bloom_bits = 16;
  opt_boundary = 3;
  opt_ceiling = 0;
  opt_differences = 1;
  opt_fastidious = 0;
  opt_gap_extension_penalty = 4;
  opt_gap_opening_penalty = 12;
  opt_help = 0;
  opt_internal_structure = 0;
  opt_log = 0;
  opt_match_reward = 5;
  opt_mismatch_penalty = 4;
  opt_mothur = 0;
  opt_no_otu_breaking = 0;
  opt_output_file = 0;
  opt_seeds = 0;
  opt_statistics_file = 0;
  opt_threads = 1;
  opt_uclust_file = 0;
  opt_usearch_abundance = 0;
  opt_version = 0;

  opterr = 1;

  char short_options[] = "a:b:c:d:e:fg:hi:l:m:no:p:rs:t:u:vw:y:z";

  /* unused short option letters: jkqx */

  static struct option long_options[] =
  {
    {"append-abundance",      required_argument, NULL, 'a' },
    {"boundary",              required_argument, NULL, 'b' },
    {"ceiling",               required_argument, NULL, 'c' },
    {"differences",           required_argument, NULL, 'd' },
    {"gap-extension-penalty", required_argument, NULL, 'e' },
    {"fastidious",            no_argument,       NULL, 'f' },
    {"gap-opening-penalty",   required_argument, NULL, 'g' },
    {"help",                  no_argument,       NULL, 'h' },
    {"internal-structure",    required_argument, NULL, 'i' },
    {"log",                   required_argument, NULL, 'l' },
    {"match-reward",          required_argument, NULL, 'm' },
    {"no-otu-breaking",       no_argument,       NULL, 'n' },
    {"output-file",           required_argument, NULL, 'o' },
    {"mismatch-penalty",      required_argument, NULL, 'p' },
    {"mothur",                no_argument,       NULL, 'r' },
    {"statistics-file",       required_argument, NULL, 's' },
    {"threads",               required_argument, NULL, 't' },
    {"uclust-file",           required_argument, NULL, 'u' },
    {"version",               no_argument,       NULL, 'v' },
    {"seeds",                 required_argument, NULL, 'w' },
    {"bloom-bits",            required_argument, NULL, 'y' },
    {"usearch-abundance",     no_argument,       NULL, 'z' },
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

    /* check if any option is specified more than once */

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
                    "Error: Option -%c or --%s specified more than once.\n",
                    c,
                    long_options[longoptindex].name);
            exit(1);
          }
        used_options[optindex] = 1;
      }

    switch(c)
      {
      case 'a':
        /* append-abundance */
        opt_append_abundance = args_long(optarg, "-a or --append-abundance");
        break;

      case 'b':
        /* boundary */
        opt_boundary = args_long(optarg, "-b or --boundary");
        break;

      case 'c':
        /* ceiling */
        opt_ceiling = args_long(optarg, "-c or --ceiling");
        break;

      case 'd':
        /* differences (resolution) */
        opt_differences = args_long(optarg, "-d or --differences");
        break;

      case 'e':
        /* gap extension penalty */
        opt_gap_extension_penalty = args_long(optarg, "-e or --gap-extension-penalty");
        break;

      case 'f':
        /* fastidious */
        opt_fastidious = 1;
        break;

      case 'g':
        /* gap-opening-penalty */
        opt_gap_opening_penalty = args_long(optarg, "-g or --gap-opening-penalty");
        break;

      case 'h':
        /* help */
        opt_help = 1;
        break;

      case 'i':
        /* internal-structure */
        opt_internal_structure = optarg;
        break;

      case 'l':
        /* log */
        opt_log = optarg;
        break;

      case 'm':
        /* match-reward */
        opt_match_reward = args_long(optarg, "-m or --match-reward");
        break;

      case 'n':
        /* no-otu-breaking */
        opt_no_otu_breaking = 1;
        break;

      case 'o':
        /* output-file */
        opt_output_file = optarg;
        break;

      case 'p':
        /* mismatch-penalty */
        opt_mismatch_penalty = args_long(optarg, "-p or --mismatch-penalty");
        break;

      case 'r':
        /* mothur */
        opt_mothur = 1;
        break;

      case 's':
        /* statistics-file */
        opt_statistics_file = optarg;
        break;

      case 't':
        /* threads */
        opt_threads = args_long(optarg, "-t or --threads");
        break;

      case 'u':
        /* uclust-file */
        opt_uclust_file = optarg;
        break;

      case 'v':
        /* version */
        opt_version = 1;
        break;

      case 'w':
        /* seeds */
        opt_seeds = optarg;
        break;

      case 'y':
        /* bloom-bits */
        opt_bloom_bits = args_long(optarg, "-y or --bloom-bits");
        break;

      case 'z':
        /* usearch-abundance */
        opt_usearch_abundance = 1;
        break;

      default:
        show_header();
        args_usage();
        exit(1);
        break;
    }
  }

  if (optind < argc)
    input_filename = argv[optind];

  if ((opt_threads < 1) || (opt_threads > MAX_THREADS))
    {
      fprintf(stderr, "\nError: Illegal number of threads specified with -t or --threads, must be in the range 1 to %d.\n", MAX_THREADS);
      exit(1);
    }

  if ((opt_differences < 0) || (opt_differences > 255))
    fatal("Illegal number of differences specified with -d or --differences, must be in the range 0 to 255.");

  if (opt_fastidious && (opt_differences != 1))
    fatal("Fastidious mode (specified with -f or --fastidious) only works when the resolution (specified with -d or --differences) is 1.");

  if (!opt_fastidious)
    {
      if (used_options[1])
        fatal("Option -b or --boundary specified without -f or --fastidious.\n");
      if (used_options[2])
        fatal("Option -c or --ceiling specified without -f or --fastidious.\n");
      if (used_options[24])
        fatal("Option -y or --bloom-bits specified without -f or --fastidious.\n");
    }

  if (opt_differences < 2)
    {
      if (used_options[12])
        fatal("Option -m or --match-reward specified when d < 2.");
      if (used_options[15])
        fatal("Option -p or --mismatch-penalty specified when d < 2.");
      if (used_options[6])
        fatal("Option -g or --gap-opening-penalty specified when d < 2.");
      if (used_options[4])
        fatal("Option -e or --gap-extension-penalty specified when d < 2.");
    }

  if (opt_gap_opening_penalty < 0)
    fatal("Illegal gap opening penalty specified with -g or --gap-opening-penalty, must not be negative.");

  if (opt_gap_extension_penalty < 0)
    fatal("Illegal gap extension penalty specified with -e or --gap-extension-penalty, must not be negative.");

  if ((opt_gap_opening_penalty + opt_gap_extension_penalty) < 1)
    fatal("Illegal gap penalties specified, the sum of the gap open and the gap extension penalty must be at least 1.");

  if (opt_match_reward < 1)
    fatal("Illegal match reward specified with -m or --match-reward, must be at least 1.");

  if (opt_mismatch_penalty < 1)
    fatal("Illegal mismatch penalty specified with -p or --mismatch-penalty, must be at least 1.");

  if (opt_boundary < 2)
    fatal("Illegal boundary specified with -b or --boundary, must be at least 2.");

  if (used_options[2] && ((opt_ceiling < 8) || (opt_ceiling > 1073741824)))
    fatal("Illegal memory ceiling specified with -c or --ceiling, must be in the range 8 to 1073741824 MB.");

  if ((opt_bloom_bits < 2) || (opt_bloom_bits > 64))
    fatal("Illegal number of Bloom filter bits specified with -y or --bloom-bits, must be in the range 2 to 64.");

  if (used_options[0] && (opt_append_abundance < 1))
    fatal("Illegal abundance value specified with -a or --append-abundance, must be at least 1.");


  /* replace filename "-" by "/dev/stdin" for input file options */

  if (!strcmp(input_filename, DASH_FILENAME))
    input_filename = STDIN_NAME;

  /* replace filename "-" by "/dev/stdout" for output file options */

  char * * stdout_options[] =
    {
      & opt_internal_structure,
      & opt_log,
      & opt_output_file,
      & opt_statistics_file,
      & opt_uclust_file,
      & opt_seeds,
      0
    };

  int o = 0;
  while(char * * stdout_opt = stdout_options[o++])
    if ((*stdout_opt) && (!strcmp(*stdout_opt, DASH_FILENAME)))
      *stdout_opt = STDOUT_NAME;
}

void open_files()
{
  /* open files */

  if (opt_log)
    {
      logfile = fopen(opt_log, "w");
      if (! logfile)
        fatal("Unable to open log file for writing.");
    }
  else
    logfile = stderr;

  if (opt_output_file)
    {
      outfile = fopen(opt_output_file, "w");
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

  if (opt_statistics_file)
    {
      statsfile = fopen(opt_statistics_file, "w");
      if (! statsfile)
        fatal("Unable to open statistics file for writing.");
    }
  else
    statsfile = 0;

  if (opt_uclust_file)
    {
      uclustfile = fopen(opt_uclust_file, "w");
      if (! uclustfile)
        fatal("Unable to open uclust file for writing.");
    }
  else
    uclustfile = 0;

  if (opt_internal_structure)
    {
      internal_structure_file = fopen(opt_internal_structure, "w");
      if (! internal_structure_file)
        fatal("Unable to open internal structure file for writing.");
    }
  else
    internal_structure_file = 0;
}

void close_files()
{
  if (opt_internal_structure)
    fclose(internal_structure_file);

  if (uclustfile)
    fclose(uclustfile);

  if (statsfile)
    fclose(statsfile);

  if (opt_seeds)
    fclose(fp_seeds);

  if (opt_output_file)
    fclose(outfile);

  if (opt_log)
    fclose(logfile);
}

int main(int argc, char** argv)
{
#ifdef __x86_64__
  cpu_features_detect();

  if (!sse2_present)
    fatal("This program requires a processor with SSE2 instructions.\n");
#endif

  srandom(1);

  args_init(argc, argv);

  open_files();

  if (opt_version || opt_help)
    {
      show_header();
      if (opt_help)
        args_usage();
      close_files();
      exit(0);
    }

  penalty_mismatch = 2 * opt_match_reward + 2 * opt_mismatch_penalty;
  penalty_gapopen = 2 * opt_gap_opening_penalty;
  penalty_gapextend = opt_match_reward + 2 * opt_gap_extension_penalty;

  penalty_factor = gcd(gcd(penalty_mismatch, penalty_gapopen), penalty_gapextend);

  penalty_mismatch /= penalty_factor;
  penalty_gapopen /= penalty_factor;
  penalty_gapextend /= penalty_factor;

  show_header();

  args_show();

  fprintf(logfile, "\n");

  db_read(input_filename);

  fprintf(logfile, "Database info:     %lu nt", db_getnucleotidecount());
  fprintf(logfile, " in %lu sequences,", db_getsequencecount());
  fprintf(logfile, " longest %lu nt\n", db_getlongestsequence());

  dbsequencecount = db_getsequencecount();

  score_matrix_init();

  switch (opt_differences)
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

  score_matrix_free();

  db_free();

  close_files();
}
