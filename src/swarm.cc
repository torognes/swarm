/*
    SWARM

    Copyright (C) 2012-2020 Torbjorn Rognes and Frederic Mahe

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

static char * progname;
char * input_filename;

char * opt_internal_structure;
char * opt_log;
char * opt_network_file;
char * opt_output_file;
char * opt_seeds;
char * opt_statistics_file;
char * opt_uclust_file;

int64_t opt_append_abundance;
int64_t opt_bloom_bits;
int64_t opt_boundary;
int64_t opt_ceiling;
bool opt_disable_sse3 {false};
int64_t opt_differences;
bool opt_fastidious {false};
int64_t opt_gap_extension_penalty;
int64_t opt_gap_opening_penalty;
bool opt_help {false};
int64_t opt_match_reward;
int64_t opt_mismatch_penalty;
bool opt_mothur {false};
int64_t opt_no_otu_breaking;
int64_t opt_threads;
int64_t opt_usearch_abundance;
int64_t opt_version;

int64_t penalty_factor;
int64_t penalty_gapextend;
int64_t penalty_gapopen;
int64_t penalty_mismatch;

/* Other variables */

int64_t mmx_present {0};
int64_t sse_present {0};
int64_t sse2_present {0};
int64_t sse3_present {0};
int64_t ssse3_present {0};
int64_t sse41_present {0};
int64_t sse42_present {0};
int64_t popcnt_present {0};
int64_t avx_present {0};
int64_t avx2_present {0};

static uint64_t dbsequencecount {0};

uint64_t duplicates_found {0};

FILE * outfile {nullptr};
FILE * statsfile {nullptr};
FILE * uclustfile {nullptr};
FILE * logfile {nullptr};
FILE * internal_structure_file {nullptr};
FILE * fp_seeds {nullptr};
FILE * network_file {nullptr};

char sym_nt[] = "-ACGT                           ";

static char dash[] = "-";
static char * DASH_FILENAME = dash;

#ifdef __x86_64__

void cpuid(unsigned int f1,
           unsigned int f2,
           unsigned int & a,
           unsigned int & b,
           unsigned int & c,
           unsigned int & d);

void cpu_features_detect();

void cpu_features_show();

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
  unsigned int a {0};
  unsigned int b {0};
  unsigned int c {0};
  unsigned int d {0};

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
  if (mmx_present){
    fprintf(logfile, " mmx");
  }
  if (sse_present) {
    fprintf(logfile, " sse");
  }
  if (sse2_present) {
    fprintf(logfile, " sse2");
  }
  if (sse3_present) {
    fprintf(logfile, " sse3");
  }
  if (ssse3_present) {
    fprintf(logfile, " ssse3");
  }
  if (sse41_present) {
    fprintf(logfile, " sse4.1");
  }
  if (sse42_present) {
    fprintf(logfile, " sse4.2");
  }
  if (popcnt_present) {
    fprintf(logfile, " popcnt");
  }
  if (avx_present) {
    fprintf(logfile, " avx");
  }
  if (avx2_present) {
    fprintf(logfile, " avx2");
  }
  fprintf(logfile, "\n");
}

#endif

int64_t args_long(char * str, const char * option);
void args_show();
void args_usage();
void show_header();
void args_init(int argc, char **argv);
void open_files();
void close_files();

int64_t args_long(char * str, const char * option)
{
  char * endptr {nullptr};
  int64_t temp = strtol(str, & endptr, 10);
  if (*endptr)
    {
      fprintf(stderr, "\nInvalid numeric argument for option %s\n", option);
      exit(1);
    }
  return temp;
}

void args_show()
{
#ifdef __x86_64__
  cpu_features_show();
#endif

  fprintf(logfile, "Database file:     %s\n", input_filename);
  fprintf(logfile, "Output file:       %s\n", opt_output_file);
  if (opt_statistics_file) {
    fprintf(logfile, "Statistics file:   %s\n", opt_statistics_file);
  }
  if (opt_uclust_file) {
    fprintf(logfile, "Uclust file:       %s\n", opt_uclust_file);
  }
  if (opt_internal_structure) {
    fprintf(logfile, "Int. struct. file  %s\n", opt_internal_structure);
  }
  if (opt_network_file) {
    fprintf(logfile, "Network file       %s\n", opt_network_file);
  }
  fprintf(logfile, "Resolution (d):    %" PRId64 "\n", opt_differences);
  fprintf(logfile, "Threads:           %" PRId64 "\n", opt_threads);

  if (opt_differences > 1)
    {
      fprintf(logfile,
              "Scores:            match: %" PRId64 ", mismatch: %" PRId64 "\n",
              opt_match_reward, opt_mismatch_penalty);
      fprintf(logfile,
              "Gap penalties:     opening: %" PRId64 ", extension: %" PRId64 "\n",
              opt_gap_opening_penalty, opt_gap_extension_penalty);
      fprintf(logfile,
              "Converted costs:   mismatch: %" PRId64 ", gap opening: %" PRId64 ", "
              "gap extension: %" PRId64 "\n",
              penalty_mismatch, penalty_gapopen, penalty_gapextend);
    }
  fprintf(logfile, "Break OTUs:        %s\n",
          opt_no_otu_breaking ? "No" : "Yes");
  if (opt_fastidious) {
    fprintf(logfile, "Fastidious:        Yes, with boundary %" PRId64 "\n",
            opt_boundary);
  }
  else {
    fprintf(logfile, "Fastidious:        No\n");
  }
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
  fprintf(stderr, " -x, --disable-sse3                  disable SSE3 and later x86 instructions\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Clustering options:\n");
  fprintf(stderr, " -d, --differences INTEGER           resolution (1)\n");
  fprintf(stderr, " -n, --no-otu-breaking               never break OTUs (not recommended!)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Fastidious options (only when d = 1):\n");
  fprintf(stderr, " -b, --boundary INTEGER              min mass of large OTUs (3)\n");
  fprintf(stderr, " -c, --ceiling INTEGER               max memory in MB for Bloom filter (unlim.)\n");
  fprintf(stderr, " -f, --fastidious                    link nearby low-abundance swarms\n");
  fprintf(stderr, " -y, --bloom-bits INTEGER            bits used per Bloom filter entry (16)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Input/output options:\n");
  fprintf(stderr, " -a, --append-abundance INTEGER      value to use when abundance is missing\n");
  fprintf(stderr, " -i, --internal-structure FILENAME   write internal OTU structure to file\n");
  fprintf(stderr, " -j, --network_file FILENAME         dump sequence network to file\n");
  fprintf(stderr, " -l, --log FILENAME                  log to file, not to stderr\n");
  fprintf(stderr, " -o, --output-file FILENAME          output result to file (stdout)\n");
  fprintf(stderr, " -r, --mothur                        output using mothur-like format\n");
  fprintf(stderr, " -s, --statistics-file FILENAME      dump OTU statistics to file\n");
  fprintf(stderr, " -u, --uclust-file FILENAME          output using UCLUST-like format to file\n");
  fprintf(stderr, " -w, --seeds FILENAME                write OTU representatives to FASTA file\n");
  fprintf(stderr, " -z, --usearch-abundance             abundance annotation in usearch style\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Pairwise alignment advanced options (only when d > 1):\n");
  fprintf(stderr, " -m, --match-reward INTEGER          reward for nucleotide match (5)\n");
  fprintf(stderr, " -p, --mismatch-penalty INTEGER      penalty for nucleotide mismatch (4)\n");
  fprintf(stderr, " -g, --gap-opening-penalty INTEGER   gap open penalty (12)\n");
  fprintf(stderr, " -e, --gap-extension-penalty INTEGER gap extension penalty (4)\n");

#ifndef __WIN32
  fprintf(stderr, "\n");
  fprintf(stderr, "See 'man swarm' for more details.\n");
#endif
}

void show_header()
{
  fprintf(logfile, "%s\n", swarm_version);
  fprintf(logfile, "Copyright (C) 2012-2020 Torbjorn Rognes and Frederic Mahe\n");
  fprintf(logfile, "https://github.com/torognes/swarm\n");
  fprintf(logfile, "\n");
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
  opt_disable_sse3 = false;
  opt_fastidious = false;
  opt_gap_extension_penalty = 4;
  opt_gap_opening_penalty = 12;
  opt_help = false;
  opt_internal_structure = nullptr;
  opt_log = nullptr;
  opt_match_reward = 5;
  opt_mismatch_penalty = 4;
  opt_mothur = false;
  opt_network_file = nullptr;
  opt_no_otu_breaking = 0;
  opt_output_file = DASH_FILENAME;
  opt_seeds = nullptr;
  opt_statistics_file = nullptr;
  opt_threads = 1;
  opt_uclust_file = nullptr;
  opt_usearch_abundance = 0;
  opt_version = 0;
  constexpr unsigned int min_ceiling {8};
  constexpr unsigned int max_ceiling {1 << 30};  // 1,073,741,824 (MiB of RAM)

  opterr = 1;

  char short_options[] = "a:b:c:d:e:fg:hi:j:l:m:no:p:rs:t:u:vw:xy:z";

  /* unused short option letters: kq */

  static struct option long_options[] =
  {
    {"append-abundance",      required_argument, nullptr, 'a' },
    {"boundary",              required_argument, nullptr, 'b' },
    {"ceiling",               required_argument, nullptr, 'c' },
    {"differences",           required_argument, nullptr, 'd' },
    {"gap-extension-penalty", required_argument, nullptr, 'e' },
    {"fastidious",            no_argument,       nullptr, 'f' },
    {"gap-opening-penalty",   required_argument, nullptr, 'g' },
    {"help",                  no_argument,       nullptr, 'h' },
    {"internal-structure",    required_argument, nullptr, 'i' },
    {"log",                   required_argument, nullptr, 'l' },
    {"network-file",          required_argument, nullptr, 'j' },
    {"match-reward",          required_argument, nullptr, 'm' },
    {"no-otu-breaking",       no_argument,       nullptr, 'n' },
    {"output-file",           required_argument, nullptr, 'o' },
    {"mismatch-penalty",      required_argument, nullptr, 'p' },
    {"mothur",                no_argument,       nullptr, 'r' },
    {"statistics-file",       required_argument, nullptr, 's' },
    {"threads",               required_argument, nullptr, 't' },
    {"uclust-file",           required_argument, nullptr, 'u' },
    {"version",               no_argument,       nullptr, 'v' },
    {"seeds",                 required_argument, nullptr, 'w' },
    {"disable-sse3",          no_argument,       nullptr, 'x' },
    {"bloom-bits",            required_argument, nullptr, 'y' },
    {"usearch-abundance",     no_argument,       nullptr, 'z' },
    {nullptr,                 0,                 nullptr, 0 }
  };

  std::array<int, 26> used_options {0};  // set int values to zero by default

  int option_index {0};
  int c {0};

  while ((c = getopt_long(argc, argv, short_options, long_options, &option_index)) != -1)
  {

    /* check if any option is specified more than once */

    if ((c >= 'a') && (c <= 'z'))
      {
        int optindex = c - 'a';
        if (used_options[optindex] == 1)
          {
            int longoptindex {0};
            while (long_options[longoptindex].name)
              {
                if (long_options[longoptindex].val == c) {
                  break;
                }
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
        opt_fastidious = true;
        break;

      case 'g':
        /* gap-opening-penalty */
        opt_gap_opening_penalty = args_long(optarg, "-g or --gap-opening-penalty");
        break;

      case 'h':
        /* help */
        opt_help = true;
        break;

      case 'i':
        /* internal-structure */
        opt_internal_structure = optarg;
        break;

      case 'j':
        /* network-file */
        opt_network_file = optarg;
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
        opt_mothur = true;
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

      case 'x':
        /* disable-sse3 */
        opt_disable_sse3 = true;
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
    }
  }

  if (optind < argc) {
    input_filename = argv[optind];
  }

  if ((opt_threads < 1) || (opt_threads > max_threads))
    {
      fprintf(stderr, "\nError: Illegal number of threads specified with "
              "-t or --threads, must be in the range 1 to %u.\n", max_threads);
      exit(1);
    }

  if ((opt_differences < 0) || (opt_differences > UINT8_MAX)) {
    fatal("Illegal number of differences specified with -d or --differences, "
          "must be in the range 0 to 255.");
  }

  if (opt_fastidious && (opt_differences != 1)) {
    fatal("Fastidious mode (specified with -f or --fastidious) only works "
          "when the resolution (specified with -d or --differences) is 1.");
  }

  if (opt_disable_sse3 && (opt_differences < 2)) {
    fatal("Option --disable-sse3 or -x has no effect when d < 2. " 
          "(SSE3 instructions are only used when d > 1).");
  }

  if (!opt_fastidious)
    {
      if (used_options[1]) {
        fatal("Option -b or --boundary specified without -f or --fastidious.\n");
      }
      if (used_options[2]) {
        fatal("Option -c or --ceiling specified without -f or --fastidious.\n");
      }
      if (used_options[24]) {
        fatal("Option -y or --bloom-bits specified without -f or --fastidious.\n");
      }
    }

  if (opt_differences < 2)
    {
      if (used_options[12]) {
        fatal("Option -m or --match-reward specified when d < 2.");
      }
      if (used_options[15]) {
        fatal("Option -p or --mismatch-penalty specified when d < 2.");
      }
      if (used_options[6]) {
        fatal("Option -g or --gap-opening-penalty specified when d < 2.");
      }
      if (used_options[4]) {
        fatal("Option -e or --gap-extension-penalty specified when d < 2.");
      }
    }

  if (opt_gap_opening_penalty < 0) {
    fatal("Illegal gap opening penalty specified with -g or "
          "--gap-opening-penalty, must not be negative.");
  }

  if (opt_gap_extension_penalty < 0) {
    fatal("Illegal gap extension penalty specified with -e or "
          "--gap-extension-penalty, must not be negative.");
  }

  if ((opt_gap_opening_penalty + opt_gap_extension_penalty) < 1) {
    fatal("Illegal gap penalties specified, the sum of the gap open and "
          "the gap extension penalty must be at least 1.");
  }

  if (opt_match_reward < 1) {
    fatal("Illegal match reward specified with -m or --match-reward, "
          "must be at least 1.");
  }

  if (opt_mismatch_penalty < 1) {
    fatal("Illegal mismatch penalty specified with -p or --mismatch-penalty, "
          "must be at least 1.");
  }

  if (opt_boundary < 2) {
    fatal("Illegal boundary specified with -b or --boundary, "
          "must be at least 2.");
  }

  if (used_options[2] && ((opt_ceiling < min_ceiling) || (opt_ceiling > max_ceiling))) {
    fatal("Illegal memory ceiling specified with -c or --ceiling, "
          "must be in the range 8 to 1,073,741,824 MB.");
  }

  if ((opt_bloom_bits < 2) || (opt_bloom_bits > 64)) {
    fatal("Illegal number of Bloom filter bits specified with -y or "
          "--bloom-bits, must be in the range 2 to 64.");
  }

  if (used_options[0] && (opt_append_abundance < 1)) {
    fatal("Illegal abundance value specified with -a or --append-abundance, "
          "must be at least 1.");
  }

  if (opt_network_file && (opt_differences != 1)) {
    fatal("A network file can only written when d=1.");
  }
}

void open_files()
{
  /* open files */

  if (opt_log)
    {
      logfile = fopen_output(opt_log);
      if (! logfile) {
        fatal("Unable to open log file for writing.");
      }
    }

  outfile = fopen_output(opt_output_file);
  if (! outfile) {
    fatal("Unable to open output file for writing.");
  }

  if (opt_seeds)
    {
      fp_seeds = fopen_output(opt_seeds);
      if (! fp_seeds) {
        fatal("Unable to open seeds file for writing.");
      }
    }

  if (opt_statistics_file)
    {
      statsfile = fopen_output(opt_statistics_file);
      if (! statsfile) {
        fatal("Unable to open statistics file for writing.");
      }
    }

  if (opt_uclust_file)
    {
      uclustfile = fopen_output(opt_uclust_file);
      if (! uclustfile) {
        fatal("Unable to open uclust file for writing.");
      }
    }

  if (opt_internal_structure)
    {
      internal_structure_file = fopen_output(opt_internal_structure);
      if (! internal_structure_file) {
        fatal("Unable to open internal structure file for writing.");
      }
    }

  if (opt_network_file)
    {
      network_file = fopen_output(opt_network_file);
      if (! network_file) {
        fatal("Unable to open network file for writing.");
      }
    }
}

void close_files()
{
  if (network_file) {
    fclose(network_file);
  }

  if (internal_structure_file) {
    fclose(internal_structure_file);
  }

  if (uclustfile) {
    fclose(uclustfile);
  }

  if (statsfile) {
    fclose(statsfile);
  }

  if (fp_seeds) {
    fclose(fp_seeds);
  }

  if (outfile) {
    fclose(outfile);
  }

  if (logfile) {
    fclose(logfile);
  }
}

int main(int argc, char** argv)
{
  logfile = stderr;

  arch_srandom(1);

  args_init(argc, argv);

#ifdef __x86_64__
  cpu_features_detect();

  if (!sse2_present) {
    fatal("This program requires a processor with SSE2 instructions.\n");
  }

  if (opt_disable_sse3)
    {
      sse3_present = 0;
      ssse3_present = 0;
      sse41_present = 0;
      sse42_present = 0;
      popcnt_present = 0;
      avx_present = 0;
      avx2_present = 0;
    }
#endif

  open_files();

  if (opt_version || opt_help)
    {
      show_header();
      if (opt_help) {
        args_usage();
      }
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

  fprintf(logfile, "Database info:     %" PRIu64 " nt", db_getnucleotidecount());
  fprintf(logfile, " in %u sequences,", db_getsequencecount());
  fprintf(logfile, " longest %u nt\n", db_getlongestsequence());

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
