/*
  SWARM

  Copyright (C) 2012-2026 Torbjorn Rognes and Frederic Mahe

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

#include "cli.hpp"
#include "swarm.hpp"
#include "utils/fatal.hpp"
#include "utils/print_view.hpp"  // fprint, fprint_integer
#include "utils/view.hpp"  // View, make_view
#include "utils/gcd.hpp"
#include "utils/open_files.hpp"
#include "arch/x86_64/cpu_features.hpp"
#include <algorithm>  // std::min(), std::transform()
#include <array>
#include <bitset>
#include <cassert>
#include <cerrno>  // errno, ERANGE
#include <cstdint>  // int64_t
#include <cstdio>  // FILE, fclose, fputs(), stderr
#include <cstdlib>  // std::exit, std::strtoll
#include <getopt.h>  // getopt_long, optarg, optind, struct option
                     // (no_argument, required_argument)
#include <cstddef>  // std::size_t
#include <iterator>  // std::next
#include <limits>
#include <string>
// <unistd.h> is not enough: it declares plain getopt(), but getopt_long()
// is a GNU extension and lives in <getopt.h> (which mingw-w64 also
// provides, so the Windows build resolves it too).


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  // Tracks which optional arguments were explicitly specified on the
  // command line. Only options whose presence influences a downstream
  // dependency or range check need a flag here; pure boolean flags
  // (e.g. --fastidious) and always-overwritten values (e.g. file paths)
  // are read directly from Parameters.
  struct UsedOptions {
    bool append_abundance {false};
    bool boundary {false};
    bool ceiling {false};
    bool gap_extension_penalty {false};
    bool gap_opening_penalty {false};
    bool match_reward {false};
    bool mismatch_penalty {false};
    bool bloom_bits {false};
  };

  constexpr char const * swarm_version {"3.1.8"};


  /* file names and command line options */

  // Single source of truth for every command-line option. long_options
  // (the C struct option array consumed by getopt_long) and short_options
  // (the colon-encoded short-letter string) are derived from this table,
  // so adding or renaming an option only requires editing one place.
  //
  // refactoring: add option -q (no-cluster-breaking)
  // (currently unused short letters: k, q)
  struct OptionSpec {
    char short_name;
    const char * long_name;
    bool needs_arg;
  };

  constexpr std::array<OptionSpec, 24> option_specs {{
      {'a', "append-abundance",      true },
      {'b', "boundary",              true },
      {'c', "ceiling",               true },
      {'d', "differences",           true },
      {'e', "gap-extension-penalty", true },
      {'f', "fastidious",            false},
      {'g', "gap-opening-penalty",   true },
      {'h', "help",                  false},
      {'i', "internal-structure",    true },
      {'j', "network-file",          true },
      {'l', "log",                   true },
      {'m', "match-reward",          true },
      {'n', "no-otu-breaking",       false},
      {'o', "output-file",           true },
      {'p', "mismatch-penalty",      true },
      {'r', "mothur",                false},
      {'s', "statistics-file",       true },
      {'t', "threads",               true },
      {'u', "uclust-file",           true },
      {'v', "version",               false},
      {'w', "seeds",                 true },
      {'x', "disable-sse3",          false},
      {'y', "bloom-bits",            true },
      {'z', "usearch-abundance",     false},
    },};


  auto build_long_options() -> std::array<struct option, option_specs.size() + 1> {
    std::array<struct option, option_specs.size() + 1> result {{}};
    std::transform(option_specs.cbegin(), option_specs.cend(), result.begin(),
                   [](OptionSpec const & spec) -> struct option {
                     struct option converted {};
                     converted.name    = spec.long_name;
                     converted.has_arg = spec.needs_arg ? required_argument : no_argument;
                     converted.flag    = nullptr;
                     converted.val     = static_cast<unsigned char>(spec.short_name);
                     return converted;
                   });
    // last slot is the {nullptr, 0, nullptr, 0} sentinel (value-initialised
    // above, and left untouched: transform writes one entry per spec, and
    // result is one longer than option_specs)
    return result;
  }


  auto build_short_options() -> std::string {
    std::string result;
    result.reserve(option_specs.size() * 2);
    for (auto const & spec : option_specs) {
      result += spec.short_name;
      if (spec.needs_arg) {
        result += ':';
      }
    }
    return result;
  }


  std::array<struct option, option_specs.size() + 1> const long_options = build_long_options();
  std::string const short_options = build_short_options();


  std::array<char const *, 18> const header_message {{
      "Swarm ", swarm_version,
      "\n",
      "Copyright (C) 2012-2026 Torbjorn Rognes and Frederic Mahe\n",
      "https://github.com/torognes/swarm\n",
      "\n",
      "Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M (2014)\n",
      "Swarm: robust and fast clustering method for amplicon-based studies\n",
      "PeerJ 2:e593 https://doi.org/10.7717/peerj.593\n",
      "\n",
      "Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M (2015)\n",
      "Swarm v2: highly-scalable and high-resolution amplicon clustering\n",
      "PeerJ 3:e1420 https://doi.org/10.7717/peerj.1420\n",
      "\n",
      "Mahe F, Czech L, Stamatakis A, Quince C, de Vargas C, Dunthorn M, Rognes T (2022)\n",
      "Swarm v3: towards tera-scale amplicon clustering\n",
      "Bioinformatics 38:1, 267-269 https://doi.org/10.1093/bioinformatics/btab493\n",
      "\n",
    },};


#ifdef _WIN32
  constexpr std::size_t args_usage_count {36};
#else
  constexpr std::size_t args_usage_count {38};
#endif

  std::array<char const *, args_usage_count> const args_usage_message {{
      /*0         1         2         3         4         5         6         7          */
      /*01234567890123456789012345678901234567890123456789012345678901234567890123456789 */
      "Usage: swarm [OPTIONS] [FASTAFILE]\n",
      "\n",
      "General options:\n",
      " -h, --help                          display this help and exit\n",
      " -t, --threads INTEGER               number of threads to use (1)\n",
      " -v, --version                       display version information and exit\n",
      "\n",
      "Clustering options:\n",
      " -d, --differences INTEGER           resolution (1)\n",
      " -n, --no-otu-breaking               never break clusters (not recommended!)\n",
      "\n",
      "Fastidious options (only when d = 1):\n",
      " -b, --boundary INTEGER              min mass of large clusters (3)\n",
      " -c, --ceiling INTEGER               max memory in MB for Bloom filter (unlim.)\n",
      " -f, --fastidious                    link nearby low-abundance swarms\n",
      " -y, --bloom-bits INTEGER            bits used per Bloom filter entry (16)\n",
      "\n",
      "Input/output options:\n",
      " -a, --append-abundance INTEGER      value to use when abundance is missing\n",
      " -i, --internal-structure FILENAME   write internal cluster structure to file\n",
      " -j, --network-file FILENAME         dump sequence network to file\n",
      " -l, --log FILENAME                  log to file, not to stderr\n",
      " -o, --output-file FILENAME          output result to file (stdout)\n",
      " -r, --mothur                        output using mothur-like format\n",
      " -s, --statistics-file FILENAME      dump cluster statistics to file\n",
      " -u, --uclust-file FILENAME          output using UCLUST-like format to file\n",
      " -w, --seeds FILENAME                write cluster representatives to FASTA file\n",
      " -z, --usearch-abundance             abundance annotation in usearch style\n",
      "\n",
      "Pairwise alignment advanced options (only when d > 1):\n",
      " -m, --match-reward INTEGER          reward for nucleotide match (5)\n",
      " -p, --mismatch-penalty INTEGER      penalty for nucleotide mismatch (4)\n",
      " -g, --gap-opening-penalty INTEGER   gap open penalty (12)\n",
      " -e, --gap-extension-penalty INTEGER gap extension penalty (4)\n",
      " -x, --disable-sse3                  disable SSE3 and later x86 instructions\n",
#ifndef _WIN32
      "\n",
      "See 'man swarm' for more details.\n",
#endif
      "\n",
    },};


  // std::strtoll rather than std::stoll: stoll reports failure by throwing,
  // and swarm handles no exceptions. std::from_chars is C++17. So the C
  // function with errno checked around it is the only correct option at
  // this standard level, not a leftover -- hence the errno = 0 before the
  // call and the three conditions after it.
  // C++17 refactoring: std::from_chars
  auto args_long(char const * const str, char const * const option) -> int64_t {
    static constexpr int base_value {10};
    char * endptr {nullptr};
    errno = 0;
    auto const number = std::strtoll(str, &endptr, base_value);
    bool const empty_input {endptr == str};
    bool const trailing_garbage {*endptr != '\0'};
    bool const out_of_range {errno == ERANGE};
    if (empty_input or trailing_garbage or out_of_range)
      {
        fatal("Invalid numeric argument for option ", option, ".\n\n",
              "Frequent causes are:\n",
              " - a missing space between an argument and the next option,\n",
              " - a long option name not starting with a double dash\n",
              "   (swarm accepts '--help' or '-h', but not '-help')\n\n",
              "Please see 'swarm --help' for more details.");
      }
    return static_cast<int64_t>(number);
  }


  template <std::size_t N>
  auto show(std::array<char const *, N> const & message,
            std::FILE * const log_stream) -> void {
    for (char const * message_element : message) {
      static_cast<void>(std::fputs(message_element, log_stream));
    }
  }

  auto show_header_message(std::FILE * const log_stream) -> void {
    show(header_message, log_stream);
  }


  auto show_help_or_version_and_exit(struct Parameters const & parameters) -> void {
    if (parameters.opt_version) {
      show(header_message, parameters.logfile);
      std::exit(EXIT_SUCCESS);
    }
    if (parameters.opt_help) {
      show(header_message, parameters.logfile);
      show(args_usage_message, parameters.logfile);
      std::exit(EXIT_SUCCESS);
    }
  }


  auto args_show(struct Parameters const & parameters) -> void {
#ifdef __x86_64__
    cpu_features_show(parameters);
#endif

    // Every line here is a padded label followed by one value. The file
    // names arrive as std::string, so they go out as a View over their own
    // bytes rather than through a "%s" and a c_str().
    auto * const log = parameters.logfile;

    fprint(log, "Database file:     ");
    fprint(log, make_view(parameters.input_filename));
    fprint(log, '\n');

    fprint(log, "Output file:       ");
    fprint(log, make_view(parameters.opt_output_file));
    fprint(log, '\n');

    if (not parameters.opt_statistics_file.empty()) {
      fprint(log, "Statistics file:   ");
      fprint(log, make_view(parameters.opt_statistics_file));
      fprint(log, '\n');
    }
    if (not parameters.opt_uclust_file.empty()) {
      fprint(log, "Uclust file:       ");
      fprint(log, make_view(parameters.opt_uclust_file));
      fprint(log, '\n');
    }
    if (not parameters.opt_internal_structure.empty()) {
      fprint(log, "Int. struct. file: ");
      fprint(log, make_view(parameters.opt_internal_structure));
      fprint(log, '\n');
    }
    if (not parameters.opt_network_file.empty()) {
      fprint(log, "Network file:      ");
      fprint(log, make_view(parameters.opt_network_file));
      fprint(log, '\n');
    }

    fprint(log, "Resolution (d):    ");
    fprint_integer(log, parameters.opt_differences);
    fprint(log, '\n');

    fprint(log, "Threads:           ");
    fprint_integer(log, parameters.opt_threads.count());
    fprint(log, '\n');

    if (parameters.opt_differences > 1)
      {
        fprint(log, "Scores:            match: ");
        fprint_integer(log, parameters.opt_match_reward);
        fprint(log, ", mismatch: ");
        fprint_integer(log, parameters.opt_mismatch_penalty);
        fprint(log, '\n');

        fprint(log, "Gap penalties:     opening: ");
        fprint_integer(log, parameters.opt_gap_opening_penalty);
        fprint(log, ", extension: ");
        fprint_integer(log, parameters.opt_gap_extension_penalty);
        fprint(log, '\n');

        fprint(log, "Converted costs:   mismatch: ");
        fprint_integer(log, parameters.penalty_mismatch);
        fprint(log, ", gap opening: ");
        fprint_integer(log, parameters.penalty_gapopen);
        fprint(log, ", gap extension: ");
        fprint_integer(log, parameters.penalty_gapextend);
        fprint(log, '\n');
      }

    // Two branches rather than a ternary: 'cond ? "No\n" : "Yes\n"' has
    // operands of different array types, so it decays to a char const * and
    // matches none of the fprint overloads.
    fprint(log, "Break clusters:    ");
    if (parameters.opt_no_cluster_breaking) {
      fprint(log, "No\n");
    }
    else {
      fprint(log, "Yes\n");
    }

    if (parameters.opt_fastidious) {
      fprint(log, "Fastidious:        Yes, with boundary ");
      fprint_integer(log, parameters.opt_boundary);
      fprint(log, '\n');
    }
    else {
      fprint(log, "Fastidious:        No\n");
    }
    fprint(log, '\n');
  }


  auto fatal_duplicate_option(int const option_character) -> void {
    // Find the matching long option name to include in the error message.
    char const * long_name = "";
    for (auto const & long_option : long_options) {
      // the sentinel entry is never reached: option_character always
      // matches a short option derived from option_specs
      assert(long_option.name != nullptr);
      if (long_option.val == option_character) {
        long_name = long_option.name;
        break;
      }
    }
    fatal("Option -", static_cast<char>(option_character),
          " or --", long_name, " specified more than once.");
  }


  // Per-option validators: range-check the signed value returned by
  // args_long() and narrow it to the type the option is stored as. They
  // run at parse time, before the value is narrowed, so out-of-range or
  // negative input is rejected before it can wrap around. Cross-option
  // checks (e.g. option combinations) stay in the validate_* functions
  // called later from args_check().
  auto validate_threading(int64_t const threads) -> ThreadCount {
    // copied into locals of the argument's own type: ThreadCount owns the
    // range, but fatal() binds its arguments to references, and a static
    // constexpr data member has no out-of-line definition to bind to in C++11
    static constexpr int64_t min_threads {ThreadCount::minimum};
    static constexpr int64_t max_threads {ThreadCount::maximum};
    if ((threads < min_threads) or (threads > max_threads)) {
      fatal("Illegal number of threads specified with "
            "-t or --threads, must be in the range ", min_threads,
            " to ", max_threads, ".");
    }
    return ThreadCount{static_cast<std::uint32_t>(threads)};
  }


  auto validate_clustering(int64_t const differences) -> std::uint64_t {
    static constexpr auto uint8_max = std::numeric_limits<uint8_t>::max();
    if ((differences < 0) or (differences > uint8_max)) {
      fatal("Illegal number of differences specified with -d or --differences, "
            "must be in the range 0 to ", static_cast<unsigned int>(uint8_max), ".");
    }
    return static_cast<std::uint64_t>(differences);
  }


  auto validate_boundary(int64_t const boundary) -> std::uint64_t {
    if (boundary < 2) {
      fatal("Illegal boundary specified with -b or --boundary, "
            "must be at least 2.");
    }
    return static_cast<std::uint64_t>(boundary);
  }


  auto validate_ceiling(int64_t const ceiling) -> std::uint64_t {
    static constexpr int64_t min_ceiling {40};
    static constexpr int64_t max_ceiling {static_cast<int64_t>(std::uint64_t{1} << 30)};  // 1,073,741,824 (MiB of RAM)
    if ((ceiling < min_ceiling) or (ceiling > max_ceiling)) {
      fatal("Illegal memory ceiling specified with -c or --ceiling, "
            "must be in the range ", min_ceiling, " to ", max_ceiling, " MB.");
    }
    return static_cast<std::uint64_t>(ceiling);
  }


  auto validate_bloom_bits(int64_t const bloom_bits) -> std::uint64_t {
    static constexpr int64_t min_bits_per_entry {2};
    static constexpr int64_t max_bits_per_entry {64};
    if ((bloom_bits < min_bits_per_entry) or (bloom_bits > max_bits_per_entry)) {
      fatal("Illegal number of Bloom filter bits specified with -y or "
            "--bloom-bits, must be in the range ", min_bits_per_entry,
            " to ", max_bits_per_entry, ".");
    }
    return static_cast<std::uint64_t>(bloom_bits);
  }


  auto args_init(int const argc, char * const * const argv, struct Parameters & parameters) -> UsedOptions {
    static constexpr std::size_t alphabet_size {26};
    UsedOptions used_options {};
    std::bitset<alphabet_size> seen_options;  // duplicate detection keyed by short letter

    while (true) {
      int option_index {0};
      int const option_character {getopt_long(argc, argv, short_options.c_str(), long_options.data(), &option_index)};

      if (option_character == -1) {
        break;
      }

      /* check if any option is specified more than once */
      if ((option_character >= 'a') and (option_character <= 'z'))
        {
          auto const bit = static_cast<std::size_t>(option_character - 'a');
          if (seen_options.test(bit)) {
            fatal_duplicate_option(option_character);
          }
          seen_options.set(bit);
        }

      switch (option_character) {
      case 'a':
        /* append-abundance */
        used_options.append_abundance = true;
        parameters.opt_append_abundance = args_long(optarg, "-a or --append-abundance");
        break;

      case 'b':
        /* boundary */
        used_options.boundary = true;
        parameters.opt_boundary = validate_boundary(args_long(optarg, "-b or --boundary"));
        break;

      case 'c':
        /* ceiling */
        used_options.ceiling = true;
        parameters.opt_ceiling = validate_ceiling(args_long(optarg, "-c or --ceiling"));
        break;

      case 'd':
        /* differences (resolution) */
        parameters.opt_differences = validate_clustering(args_long(optarg, "-d or --differences"));
        break;

      case 'e':
        /* gap extension penalty */
        used_options.gap_extension_penalty = true;
        parameters.opt_gap_extension_penalty = args_long(optarg, "-e or --gap-extension-penalty");
        break;

      case 'f':
        /* fastidious */
        parameters.opt_fastidious = true;
        break;

      case 'g':
        /* gap-opening-penalty */
        used_options.gap_opening_penalty = true;
        parameters.opt_gap_opening_penalty = args_long(optarg, "-g or --gap-opening-penalty");
        break;

      case 'h':
        /* help */
        parameters.opt_help = true;
        break;

      case 'i':
        /* internal-structure */
        parameters.opt_internal_structure = optarg;
        break;

      case 'j':
        /* network-file */
        parameters.opt_network_file = optarg;
        break;

      case 'l':
        /* log */
        parameters.opt_log = optarg;
        break;

      case 'm':
        /* match-reward */
        used_options.match_reward = true;
        parameters.opt_match_reward = args_long(optarg, "-m or --match-reward");
        break;

      case 'n':
        /* no-cluster-breaking */
        parameters.opt_no_cluster_breaking = true;
        break;

      case 'o':
        /* output-file */
        parameters.opt_output_file = optarg;
        break;

      case 'p':
        /* mismatch-penalty */
        used_options.mismatch_penalty = true;
        parameters.opt_mismatch_penalty = args_long(optarg, "-p or --mismatch-penalty");
        break;

      case 'r':
        /* mothur */
        parameters.opt_mothur = true;
        break;

      case 's':
        /* statistics-file */
        parameters.opt_statistics_file = optarg;
        break;

      case 't':
        /* threads */
        parameters.opt_threads = validate_threading(args_long(optarg, "-t or --threads"));
        break;

      case 'u':
        /* uclust-file */
        parameters.opt_uclust_file = optarg;
        break;

      case 'v':
        /* version */
        parameters.opt_version = true;
        break;

      case 'w':
        /* seeds */
        parameters.opt_seeds = optarg;
        break;

      case 'x':
        /* disable-sse3 */
        parameters.opt_disable_sse3 = true;
        break;

      case 'y':
        /* bloom-bits */
        used_options.bloom_bits = true;
        parameters.opt_bloom_bits = validate_bloom_bits(args_long(optarg, "-y or --bloom-bits"));
        break;

      case 'z':
        /* usearch-abundance */
        parameters.opt_usearch_abundance = true;
        break;

      default:
        show(header_message, parameters.logfile);
        show(args_usage_message, parameters.logfile);
        fatal();
      }
    }

    if (optind < argc) {  // external variable defined in unistd.h for
      // use with the getopt function
      parameters.input_filename = *std::next(argv, optind);
      // An empty positional argument would otherwise reach parse_fasta()
      // as an empty filename, where the "filename is always set" contract
      // is asserted. Rejected here so the contract holds by construction
      // rather than aborting a debug build on ordinary bad input.
      if (parameters.input_filename.empty()) {
        fatal("Empty input file name.");
      }
    }

#ifdef __x86_64__
    cpu_features_detect(parameters);
    cpu_features_test(parameters);
#endif

    return used_options;
  }


  auto set_alignment_scoring_system(struct Parameters & parameters) -> void {
    // Reject scoring parameters whose magnitude is large enough to
    // overflow the conversion arithmetic below (signed overflow is
    // undefined behaviour). The widest expression is (2 * reward) +
    // (2 * penalty), so any operand within +/- (INT64_MAX / 4) is safe.
    // This guard runs before args_check(): legitimate scoring values must
    // fit an 8-bit score and are orders of magnitude smaller, so this only
    // catches nonsensical input (e.g. values close to INT64_MAX) before it
    // reaches the multiplication. Out-of-range-but-non-overflowing values
    // are still rejected later by check_scoring_saturation().
    static constexpr int64_t max_scoring_value {std::numeric_limits<int64_t>::max() / 4};
    struct ScoringValue {
      int64_t value;
      char const * option;
    };
    std::array<ScoringValue, 4> const scoring_values {{
        {parameters.opt_match_reward,          "-m or --match-reward"},
        {parameters.opt_mismatch_penalty,      "-p or --mismatch-penalty"},
        {parameters.opt_gap_opening_penalty,   "-g or --gap-opening-penalty"},
        {parameters.opt_gap_extension_penalty, "-e or --gap-extension-penalty"},
      },};
    for (auto const & scoring_value : scoring_values) {
      if ((scoring_value.value > max_scoring_value) or
          (scoring_value.value < -max_scoring_value)) {
        fatal("Scoring parameter specified with ", scoring_value.option,
              " is too large.");
      }
    }

    parameters.penalty_mismatch = (2 * parameters.opt_match_reward) + (2 * parameters.opt_mismatch_penalty);
    parameters.penalty_gapopen = 2 * parameters.opt_gap_opening_penalty;
    parameters.penalty_gapextend = parameters.opt_match_reward + (2 * parameters.opt_gap_extension_penalty);

    int64_t const penalty_factor {gcd(gcd(parameters.penalty_mismatch, parameters.penalty_gapopen), parameters.penalty_gapextend)};

    // penalty_factor cannot be zero: validate_alignment() runs before this
    // function and guarantees a match reward and a mismatch penalty of at
    // least 1, so penalty_mismatch is at least 4 and the gcd is at least 1.
    parameters.penalty_mismatch /= penalty_factor;
    parameters.penalty_gapopen /= penalty_factor;
    parameters.penalty_gapextend /= penalty_factor;
  }


  auto validate_fastidious(UsedOptions const & used_options,
                           struct Parameters const & parameters) -> void {
    // Range checks for boundary, ceiling and bloom-bits run at parse time
    // (validate_boundary/validate_ceiling/validate_bloom_bits). Only the
    // cross-option checks remain here.
    if (parameters.opt_fastidious and (parameters.opt_differences != 1)) {
      fatal("Fastidious mode (specified with -f or --fastidious) only works "
            "when the resolution (specified with -d or --differences) is 1.");
    }

    if (not parameters.opt_fastidious) {
      if (used_options.boundary) {
        fatal("Option -b or --boundary specified without -f or --fastidious.");
      }
      if (used_options.ceiling) {
        fatal("Option -c or --ceiling specified without -f or --fastidious.");
      }
      if (used_options.bloom_bits) {
        fatal("Option -y or --bloom-bits specified without -f or --fastidious.");
      }
    }
  }


  auto validate_alignment(UsedOptions const & used_options,
                          struct Parameters const & parameters) -> void {
    if (parameters.opt_disable_sse3 and (parameters.opt_differences < 2)) {
      fatal("Option --disable-sse3 or -x has no effect when d < 2 "
            "(SSE3 instructions are only used when d > 1).");
    }

    if (parameters.opt_differences < 2) {
      if (used_options.match_reward) {
        fatal("Option -m or --match-reward specified when d < 2.");
      }
      if (used_options.mismatch_penalty) {
        fatal("Option -p or --mismatch-penalty specified when d < 2.");
      }
      if (used_options.gap_opening_penalty) {
        fatal("Option -g or --gap-opening-penalty specified when d < 2.");
      }
      if (used_options.gap_extension_penalty) {
        fatal("Option -e or --gap-extension-penalty specified when d < 2.");
      }
    }

    if (parameters.opt_gap_opening_penalty < 0) {
      fatal("Illegal gap opening penalty specified with -g or "
            "--gap-opening-penalty, must not be negative.");
    }

    if (parameters.opt_gap_extension_penalty < 0) {
      fatal("Illegal gap extension penalty specified with -e or "
            "--gap-extension-penalty, must not be negative.");
    }

    if ((parameters.opt_gap_opening_penalty + parameters.opt_gap_extension_penalty) < 1) {
      fatal("Illegal gap penalties specified, the sum of the gap open and "
            "the gap extension penalty must be at least 1.");
    }

    if (parameters.opt_match_reward < 1) {
      fatal("Illegal match reward specified with -m or --match-reward, "
            "must be at least 1.");
    }

    if (parameters.opt_mismatch_penalty < 1) {
      fatal("Illegal mismatch penalty specified with -p or --mismatch-penalty, "
            "must be at least 1.");
    }
  }


  auto validate_io(UsedOptions const & used_options,
                   struct Parameters const & parameters) -> void {
    if (used_options.append_abundance and (parameters.opt_append_abundance < 1)) {
      fatal("Illegal abundance value specified with -a or --append-abundance, "
            "must be at least 1.");
    }

    if ((not parameters.opt_network_file.empty()) and (parameters.opt_differences != 1)) {
      fatal("A network file can only written when d = 1.");
    }
  }


  auto check_scoring_saturation(struct Parameters const & parameters) -> void {
    static constexpr auto uint8_max = std::numeric_limits<uint8_t>::max();
    static constexpr auto uint16_max = std::numeric_limits<uint16_t>::max();
    int64_t const diff_saturation_16 = std::min((uint16_max / parameters.penalty_mismatch),
                                                (uint16_max - parameters.penalty_gapopen)
                                                / parameters.penalty_gapextend);

    // diff_saturation_16 may be negative (its uint16_max - penalty_gapopen
    // term underflows for large gap penalties), which means the scoring
    // system is already saturated; compare as signed so opt_differences
    // (0..255) correctly exceeds any negative saturation value.
    if (static_cast<int64_t>(parameters.opt_differences) > diff_saturation_16) {
      fatal("Resolution (d) too high for the given scoring system.");
    }

    if (parameters.penalty_mismatch > uint8_max) {
      fatal("Alignment scoring system yielded a mismatch penalty greater than 255, "
            "please use different parameter values.");
    }
  }


  auto args_check(UsedOptions const & used_options,
                  struct Parameters const & parameters) -> void {
    validate_fastidious(used_options, parameters);
    validate_io(used_options, parameters);
    check_scoring_saturation(parameters);
  }

}  // end of anonymous namespace


auto parse_command_line(int const argc, char * const * const argv) -> Parameters {
  Parameters parameters;
  auto const used_options = args_init(argc, argv, parameters);
  show_help_or_version_and_exit(parameters);
  validate_alignment(used_options, parameters);
  set_alignment_scoring_system(parameters);
  args_check(used_options, parameters);
  open_files(parameters);
  show_header_message(parameters.logfile);
  args_show(parameters);
  return parameters;
}
