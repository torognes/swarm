# swarm — refactoring / bug-fix plan

Derived from the critical code review of the source tree (branch
`tmp_20260626151239`). Each point below lists **selectable
implementation choices**. Pick one option per point (or amend it); the
fixes themselves are intentionally **not yet implemented** — this
document plus the matching regression tests in `swarm-tests`
(`scripts/pending_fixes.sh`, same branch name) are the deliverable for
review.

Legend for the *Test* line:
- **red test ready** — a failing regression test exists in
  `swarm-tests/scripts/pending_fixes.sh`; it turns green once the fix
  lands and should then be migrated into `fixed_bugs.sh`.
- **not black-box testable** — the defect cannot be exercised through
  the swarm CLI on this platform (build-time, foreign-OS, or
  astronomically-large input); verified by reading code / review only.

Each fix should be a small, self-contained commit with a
`Co-Authored-By: Florian Filloux <...>` trailer, checked with
`cppcheck` and the cross-compilers, per `CLAUDE.md`.

---

## BUG 1 — CLI integer divide-by-zero (SIGFPE) in scoring setup

- **Severity:** BUG (reproducible crash from the command line, exit 136).
- **Location:** `src/cli.cc:587-596` (`set_alignment_scoring_system`),
  ordering in `parse_command_line` (`src/cli.cc:723-724`).
- **Problem:** `set_alignment_scoring_system()` computes
  `penalty_factor = gcd(gcd(mismatch, gapopen), gapextend)` and divides
  the three penalties by it, **before** `args_check()` →
  `validate_alignment()` runs the positivity checks (match ≥ 1,
  mismatch ≥ 1, gaps ≥ 0, gap-sum ≥ 1). Scoring that drives all three
  penalties to zero makes `penalty_factor == 0` → integer division by
  zero. The comment at `cli.cc:593` ("would require gcd(0,0) which is
  not possible") is incorrect.
- **Repro:** `swarm -d 2 -m 0 -p 0 -e 0 -g 0 file.fa`
- **Constraint:** `check_scoring_saturation()` (inside `args_check`)
  *reads* the derived penalties, so it must stay **after**
  `set_alignment_scoring_system`. Only the `validate_alignment`
  positivity checks need to move earlier.

**Options:**

- **(1A) Reorder `validate_alignment` before scoring derivation.**
  In `parse_command_line`, call `validate_alignment` (or a split-out
  "scoring sanity" subset) before `set_alignment_scoring_system`; keep
  `validate_io`, `validate_fastidious`, and `check_scoring_saturation`
  in `args_check` after. *Pro:* removes the root cause (bad parameters
  rejected before any arithmetic), no defensive guard needed. *Con:*
  splits `args_check`, slightly reorders user-visible error precedence
  (a positivity error would now be reported before an I/O error).

- **(1B) Guard `penalty_factor == 0` inside `set_alignment_scoring_system`.**
  `if (penalty_factor == 0) { return; }` (or skip the divisions) before
  the three `/=`. *Pro:* one-line, no reordering. *Con:* masks the bad
  input rather than rejecting it; the later validator still has to fire
  to produce a correct error message; leaves the misleading comment's
  intent half-true.

- **(1C) Add explicit positivity asserts/checks at the top of
  `set_alignment_scoring_system`.** Duplicate the match≥1 / mismatch≥1 /
  gap checks (as real `fatal()` calls) immediately before computing the
  penalties. *Pro:* keeps all scoring logic in one place. *Con:*
  duplicates validation logic that already lives in
  `validate_alignment` (two places to maintain).

**Recommendation: 1A.** It eliminates the class of bug (no derived
value is ever computed from invalid parameters) without duplicating
checks. In all cases, **fix the inaccurate comment at `cli.cc:593`**.

- **Test:** red test ready — asserts the command exits with the
  validation error, not a floating-point exception (exit ≠ 136/139).

---

## BUG 2 — Fastidious: inflated "Light swarms" count + divide-by-zero with `--ceiling`

- **Severity:** BUG (always-wrong logged statistic on most fastidious
  runs; crash when combined with `--ceiling`).
- **Location:** `src/utils/algod1_statistics.cc:53-66` (counting) and
  `:118` (division); root cause `src/algod1.cc:157,290,293`.
- **Problem:** `ensure_swarm_capacity` grows `swarminfo_v` with
  `resize(size + 1024)`, so `swarminfo_v.size()` is `swarmcount` rounded
  up to a multiple of 1024, the extra entries default-constructed with
  `mass == 0`. `run_fastidious_pass` runs (`algod1.cc:290`) **before**
  the `swarminfo_v.resize(swarmcount)` at `algod1.cc:293`, so
  `count_cluster_stats` iterates the padding entries. Each padding entry
  has `mass (0) < opt_boundary`, so it is miscounted as a light swarm:
  - The logged `Light swarms: N` is overstated by the padding count
    (0–1023) whenever `swarmcount % 1024 != 0`, while the amplicon
    count beside it stays correct → inconsistent log.
  - When there are **zero genuine light swarms** but padding > 0 and at
    least one heavy swarm, the `small_clusters == 0` short-circuit
    (`algod1_fastidious.cc:547`) is defeated, so `compute_bloom_geometry`
    runs with `nucleotides_in_small_clusters == 0`; with `--ceiling`
    set, `algod1_statistics.cc:118` divides by
    `microvariants * nucleotides_in_small_clusters` → division by zero.
- **Repro (crash):** `-f --ceiling N` on input where every swarm is
  heavy (e.g. raise abundances / lower `--boundary`) and `swarmcount`
  is not a multiple of 1024.

**Options:**

- **(2A) Resize `swarminfo_v` to `swarmcount` before fastidious.**
  Move the `swarminfo_v.resize(swarmcount); shrink_to_fit();`
  (`algod1.cc:293`) to immediately after `run_clustering` /
  `overall_stats.swarmcount_adjusted = swarmcount`, i.e. before the
  `if (opt_fastidious)` block. Grafting (`attach_candidates`) works by
  index and appends no swarms, so this is safe. *Pro:* fixes both the
  count and the crash at the source; the padding never reaches any
  consumer. *Con:* a `shrink_to_fit` before fastidious may reallocate;
  negligible cost.

- **(2B) Make `count_cluster_stats` iterate only the first `swarmcount`
  entries.** It already takes `amplicon_count`; pass `swarmcount` too
  and loop `begin() … begin()+swarmcount`. *Pro:* minimal, localized to
  the statistics function; no change to vector lifetime. *Con:* leaves
  the over-sized vector around during fastidious (any *other* future
  consumer of `swarminfo_v` before the resize would hit the same trap).

- **(2C) Defensive guard only.** Add
  `if (nucleotides_in_small_clusters == 0) { /* skip bloom path */ }`
  in `run_fastidious_pass` / `compute_bloom_geometry`. *Pro:* stops the
  crash. *Con:* does **not** fix the wrong "Light swarms" count — only
  half the bug; should not be used alone.

**Recommendation: 2A** (root-cause fix for both symptoms),
optionally plus **2C** as defence-in-depth on `compute_bloom_geometry`.
2B is acceptable if you prefer to keep the resize where it is.

- **Test:** red test ready — two tests: (i) `Light swarms:` count in
  the log matches the real number of light swarms; (ii) `-f --ceiling`
  on an all-heavy dataset exits cleanly (no SIGFPE).

---

## BUG 3 — x86_64 POPCNT-fallback popcount uses MMX without EMMS + non-portable cast

- **Severity:** BUG (UB on the no-POPCNT path; flagged by two
  reviewers). Trivial fix.
- **Location:** `src/arch/x86_64/qgram_compare.cc:101`.
- **Problem:** `reinterpret_cast<uint64_t>(_mm_movepi64_pi64(vector_n))`.
  `_mm_movepi64_pi64` returns an `__m64` (MMX register): (a) the path
  enters MMX state with no following `_mm_empty()` (EMMS), so later
  x87/`double` work on the same thread can read corrupted FPU state;
  (b) `reinterpret_cast` from the `__m64` vector type to `uint64_t` is
  not a well-defined conversion. Runs on SSE2-without-POPCNT CPUs and on
  **any** machine invoked with `--disable-sse3` (which clears
  `popcnt_present`).

**Options:**

- **(3A) Use `_mm_cvtsi128_si64`.**
  `return static_cast<uint64_t>(_mm_cvtsi128_si64(vector_n));` — pure
  SSE2, extracts the low 64 bits via an integer register, no MMX, no
  EMMS. *Pro:* idiomatic, one line, removes both problems. *Con:*
  `_mm_cvtsi128_si64` is 64-bit-only — confirm no 32-bit x86 target
  (swarm targets x86_64/ARM64/PPC64/Win64, so fine).

- **(3B) Store and reload via `_mm_storel_epi64`.**
  Write the low 64 bits to a `uint64_t` with
  `_mm_storel_epi64(reinterpret_cast<__m128i *>(&result), vector_n);`.
  *Pro:* works on any SSE2 target including 32-bit. *Con:* a store/load
  round-trip; marginally slower on an already-cold path.

- **(3C) Keep MMX but add `_mm_empty()` and a proper cast.** *Pro:*
  smallest diff. *Con:* still uses MMX needlessly; EMMS is easy to
  forget on future edits — not recommended.

**Recommendation: 3A.**

- **Test:** red test ready (best-effort) — clustering result with
  `--disable-sse3` (no-POPCNT path) must equal the default-path result
  on the same input. Note: UB may not *visibly* manifest, so this test
  guards correctness of the fallback path rather than the FPU-state
  corruption directly. **not strictly black-box testable** for the EMMS
  aspect.

---

## RISK 4 — `dup()` success test rejects fd 0 and leaks it

- **Severity:** RISK (mishandled error path / fd leak; edge case).
- **Location:** `src/utils/input_output.cc:46-47` and `62-63`.
- **Problem:** `dup()` returns `-1` on failure; the code uses
  `file_descriptor > 0 ? fdopen(...) : nullptr`. If fd 0 (or 1) is
  closed, `dup()` legitimately returns 0, treated as failure and leaked
  (never `fdopen`'d / `close`'d).

**Options:**

- **(4A) Compare against the correct sentinel.** Change `> 0` to
  `>= 0` (or `!= -1`) in both spots. *Pro:* correct, minimal.
  *Con:* none.
- **(4B) Same as 4A but also name the sentinel.** Introduce
  `constexpr int invalid_fd {-1};` for readability. *Pro:* self-
  documenting, matches "avoid magic numbers". *Con:* trivially larger.

**Recommendation: 4B** (matches the project's magic-number guideline).

- **Test:** not black-box testable (requires the test harness to close
  fd 0/1 before exec; not reliably expressible as a portable bash test).

---

## RISK 5 — `unsigned int` overflow on total d=1 network edges

- **Severity:** RISK (wrong clusters / OOB on very large d=1 datasets).
- **Location:** `src/utils/algod1_internal.h:59,92`; written
  `src/utils/algod1_network.cc:144,154`; indexed `src/algod1.cc:92`,
  `src/utils/algod1_output.cc:53-55`.
- **Problem:** `Network_state::count` accumulates the **total** number
  of network edges across all amplicons and `ampinfo_s::link_start`
  indexes `network_v` (a `std::vector`), both `unsigned int`. On a
  large enough d=1 dataset the running total wraps past `UINT_MAX`,
  corrupting `link_start` offsets.

**Options:**

- **(5A) Widen to `uint64_t`.** Change `count`, `link_start`, the
  `write_network_file` count parameter, and the matching progress total
  to `uint64_t` (the natural index type for `network_v`). *Pro:*
  removes the ceiling entirely; matches the existing `uint64_t` mass
  accumulators. *Con:* slightly larger `ampinfo_s` (4 bytes/amplicon)
  — measurable memory at scale; benchmark on the d=1 paths per
  `CLAUDE.md` to confirm no regression.

- **(5B) Keep `unsigned int` but add a saturating assert/fatal.**
  Detect `count` approaching `UINT_MAX` and `fatal()` with a clear
  message. *Pro:* no memory growth. *Con:* turns a silent corruption
  into a hard stop — better than corruption, but refuses inputs that
  would otherwise work with 5A.

**Recommendation: 5A** if the per-amplicon memory cost is acceptable
(measure first); otherwise 5B as a stopgap.

- **Test:** not black-box testable (would require billions of edges).

---

## RISK 6 — 8-bit SIMD gap-boundary seed truncates instead of saturating

- **Severity:** RISK (wrong difference count with non-default,
  large gap penalties). **Touches bit-mode policy — needs a human call.**
- **Location:** `src/search8.cc:708` (and analogously
  `src/search16.cc:458`).
- **Problem:** the F0/H0 boundary seed is stored with a truncating cast
  `static_cast<BYTE>(2*(gapopen+gapextend))`. `set_bit_mode` only
  guarantees `d*(gapopen+gapextend) <= 255`, which permits
  `gapopen+gapextend > 127` at `d=1`; then `2*(go+ge)` wraps mod 256 to
  a too-small boundary cost → a gap-heavy alignment can score lower than
  it should → wrong difference count / cluster edge. The documenting
  `assert`s are compiled out under `NDEBUG`.

**Options:**

- **(6A) Tighten `set_bit_mode` to require `2*(gapopen+gapextend) <=
  255` before selecting the 8-bit kernel.** Cases that would overflow
  fall through to the 16-bit kernel. *Pro:* preserves exact scores;
  no clamping of values. *Con:* shifts a few extreme-parameter inputs to
  the slower 16-bit path; need to verify the 16-bit seed bound likewise
  holds (its analogous store at `search16.cc:458` uses a 16-bit `WORD`,
  which fits `2*(go+ge)` for all CLI-permitted values).

- **(6B) Clamp the seed: `std::min<unsigned>(255, 2*(go+ge))`.**
  *Pro:* no path change, no crash. *Con:* a clamped boundary seed is
  itself a slightly wrong score; only correct if the saturating kernel
  treats a 255-seed as "≥ threshold, skip" in all cases — must be
  verified, otherwise this trades one wrong score for another.

- **(6C) Document-only / no change.** If product decision is that such
  gap penalties are out of supported range, replace the `assert` with a
  `fatal()` in `set_bit_mode` rejecting the combination. *Pro:* honest
  about supported range. *Con:* removes functionality some users may
  rely on.

**Recommendation: 6A** (correct scores, only a performance shift for
exotic parameters) — but please confirm the desired behaviour for
`gapopen+gapextend > 127`, since this is a scoring-semantics decision.

- **Test:** red test ready (best-effort) — run d=2 with a large
  gap-open/extend so `2*(go+ge) > 255` and compare the cluster output to
  a reference computed by the slow/16-bit path; flagged for human review
  of the expected result per `CLAUDE.md`.

---

## RISK 7 — Windows memory queries ignore failure → garbage `--ceiling` sizing

- **Severity:** RISK (wrong result on Windows; the POSIX backends
  `fatal()` on failure, Windows diverges).
- **Location:** `src/os/windows/system_memory.cc:32-43`.
- **Problem:** `GetProcessMemoryInfo` / `GlobalMemoryStatusEx` return
  values are discarded and the structs are not zero-initialised; on
  failure the function returns garbage feeding Bloom-filter sizing.

**Options:**

- **(7A) Check both BOOLs and `fatal()` on failure; zero-init the
  structs.** Mirror the POSIX `fatal("Cannot determine amount of
  RAM.")`. *Pro:* parity across OSes. *Con:* none.
- **(7B) Check + fall back to a conservative default instead of
  `fatal()`.** *Pro:* never aborts. *Con:* inconsistent with the POSIX
  backends, which abort — divergent behaviour again.

**Recommendation: 7A** (cross-OS parity).

- **Test:** not black-box testable on this (Linux) platform.

---

## RISK 8 — SSSE3 shuffle guarded by `#ifdef __SSE3__` instead of `__SSSE3__`

- **Severity:** RISK (latent build/dispatch footgun; correct today only
  because of Makefile flags).
- **Location:** `src/arch/x86_64/search_dispatch.cc:47,63`,
  `src/arch/x86_64/ssse3.h:24`.
- **Problem:** `_mm_shuffle_epi8` is an SSSE3 (PSHUFB) instruction but
  is guarded by `__SSE3__`. Works because the Makefile compiles these
  TUs with `-mssse3` (which defines both macros); a `-msse3`-only build
  would mis-compile / mis-guard.

**Options:**

- **(8A) Change the guards to `#ifdef __SSSE3__`.** *Pro:* the guard
  now names the instruction set it actually protects; no behaviour
  change with current flags. *Con:* none — verify all three sites and
  re-run the cross-compilers.
- **(8B) Leave the guard, add a `static_assert`/comment** documenting
  the dependence on `-mssse3`. *Pro:* zero risk of a guard typo causing
  a different bug. *Con:* keeps the misleading macro name.

**Recommendation: 8A.**

- **Test:** not black-box testable (compile-time guard).

---

## RISK 9 — Hash-table / header sizing: latent overflow & degenerate inputs

Three independent, low-likelihood arithmetic issues, groupable into one
"sizing hardening" commit.

- **9a — `compute_hashtable_size` uint64 overflow before the assert.**
  `src/utils/hashtable_size.cc:41`: `denominator * (sequence_count + 1)`
  (denominator = 10) overflows uint64 at `sequence_count ≈ 1.84e18`,
  *below* the assert bound (~6.45e18).
- **9b — `compute_hashtable_size(0)` returns 1**, not the documented
  ≥ 2; `mask == 0` → infinite insert loop if ever reached.
- **9c — header length narrowed before the size check.**
  `src/db.cc:224-230`: `strcspn` (returns `size_t`) is narrowed to
  `unsigned int` *before* the `max_header_length` comparison, so a
  > 4 GiB header would pass the guard truncated.

**Options:**

- **(9A) Fix all three (recommended grouping).**
  - 9a: compute the 10/7 scaling in floating point, or reorder as
    `std::log(sequence_count + 1) + std::log(10.0/7.0)`, or tighten the
    assert to `< UINT64_MAX/10 - 1`.
  - 9b: `return std::max<uint64_t>(2, …);`.
  - 9c: keep the `strcspn` result as `size_t`, compare against
    `max_header_length` first, narrow only after the guard passes.
  *Pro:* removes three fragile-invariant violations cheaply.
  *Con:* none of practical consequence.

- **(9B) Fix only 9c (the crafted-input one), defer 9a/9b.** *Pro:*
  9c is the only one reachable by a (very large) input file; 9a/9b need
  ~quintillions of sequences. *Con:* leaves documented invariants
  violated.

- **(9C) Document-only for 9a/9b** (comment that the limits are far
  beyond addressable memory). *Pro:* zero code churn. *Con:* a future
  refactor could rely on the false "≥ 2" / "no overflow" invariant.

**Recommendation: 9A** (all three are tiny and remove latent traps).

- **Test:** 9c — red test ready (best-effort): a multi-MB (not 4 GiB)
  over-long header must be rejected with the header-too-long error, not
  silently accepted; this guards the *ordering* even though the wrap
  itself needs > 4 GiB. 9a/9b not black-box testable.

---

## RISK 10 — Oversized abundance reported as "missing" instead of "too large"

- **Severity:** RISK (confusing/wrong diagnostic on a pathological but
  real header).
- **Location:** `src/db.cc:461-499`.
- **Problem:** a 20-digit abundance exceeding `int64_t` makes
  `parse_abundance_digits` return `false` (ERANGE), which is
  indistinguishable from "no annotation present" → the user sees
  "Abundance annotations not found" instead of an overflow error.

**Options:**

- **(10A) Distinguish "matched digits but overflowed" from "no
  match".** Return a small status enum (or out-param) from
  `parse_abundance_digits` and emit an "abundance too large" `fatal()`
  on overflow. *Pro:* precise diagnostic. *Con:* touches the
  abundance-parsing signature (per the `feedback_signature_changes`
  memory, confirm before rippling to callers).

- **(10B) Pre-check digit count / value range and `fatal()` early.**
  When the matched digit run is ≥ 19 digits, attempt the parse and
  `fatal()` on ERANGE specifically at the call site. *Pro:* no signature
  change. *Con:* duplicates ERANGE handling at each call site
  (`find_swarm_abundance`, `find_usearch_abundance`).

**Recommendation: 10A** if you are comfortable with the
internal-signature change (it is internal, not public API); otherwise
10B. Confirm before editing the signature.

- **Test:** red test ready — a header with a 20-digit `;size=` must
  produce an overflow-specific message, not "Abundance annotations not
  found".

---

## RISK 11 — Defensive: empty / zero-length alignment paths

- **Severity:** RISK (UB in release if reachable; likely guarded
  upstream).
- **Location:** `src/utils/cigar.cc:31,60` (`input.back()` on an empty
  vector), `src/utils/nw_aligner.cc:273` (`percent_id` → `0/0` → NaN).
- **Problem:** if `backtrack()` ever returns an empty alignment, release
  builds invoke UB / produce NaN. Not confirmed reachable (empty input
  sequences appear to be rejected during parsing).

**Options:**

- **(11A) Confirm the upstream guard, add a comment, no code change.**
  *Pro:* no churn if the invariant truly holds. *Con:* relies on a
  distant invariant.
- **(11B) Add a local early-return** for `length == 0` (empty/`"0M"`
  result, `percent_id = 100.0` or 0). *Pro:* defence-in-depth, self-
  contained. *Con:* adds a branch to a hot path (negligible).

**Recommendation: 11A first** (verify the parser rejects empty
sequences — `db.cc:300`); add **11B** only if the invariant is not
firmly guaranteed.

- **Test:** red/clarifying test ready (best-effort): feed a record with
  an empty sequence and assert a clean rejection (no crash/NaN). May
  already pass on current code (documents the invariant).

---

## ENHANCEMENTS (no behaviour change)

Group into one or two cleanup commits; none alter results.

- **E1 — `-c/--ceiling` message contradicts the code.**
  `src/cli.cc:364-368` enforces a minimum of 40 but the message says
  "range 8 to …". *Options:* (E1A) fix the message to interpolate
  `min_ceiling`/`max_ceiling`; (E1B) lower the enforced minimum back to
  8 to match the message *and* the older man page — **this is a
  behaviour/man-page divergence, needs human review** per `CLAUDE.md`.
  *Recommendation:* E1A (message wrong, code matches the current man
  page). **Test: red test ready** — `swarm -c 20` error text must state
  the real lower bound (40).

- **E2 — `db.cc:917` function-local `static std::vector<char> buffer`
  in `fprintseq`.** Shared mutable state; safe only because output is
  single-threaded and one `Data` per run. *Options:* (E2A) make it a
  non-static local; (E2B) make it a per-`Data` member.
  *Recommendation:* E2A (simplest; cost negligible vs. I/O).
  *Test:* not black-box testable.

- **E3 — AVX/AVX2 reported without OSXSAVE/XGETBV check.**
  `src/arch/x86_64/cpu_features.cc:71,80`. Cosmetic today (AVX never
  gates a kernel). *Options:* (E3A) add the OSXSAVE+XGETBV check;
  (E3B) add a comment that the flags are display-only. *Recommendation:*
  E3B now, E3A only if an AVX kernel is ever added. *Test:* not
  black-box testable.

- **E4 — Minor hygiene** (single small commit): signedness of the
  gap-penalty asserts (`search8.cc:771`, `search16.cc:521`) using signed
  `numeric_limits<char/short>::max()` for unsigned `BYTE/WORD`;
  `db.cc:141` stray `;;`; `db.cc:180` misleading parameter name
  `mapped_minus_one`; `os/macos/system_memory.cc:41` C array
  `int mib[]` → `std::array`. *Recommendation:* fix all; pure cleanup.
  *Test:* not black-box testable.

---

## Suggested commit / merge order

1. **BUG 3** (3A) — trivial, isolated, removes UB.
2. **BUG 1** (1A) — crash fix + comment correction.
3. **BUG 2** (2A [+2C]) — crash + wrong-statistic fix.
4. **RISK 10** (10A/10B), **E1** (E1A) — user-facing diagnostics.
5. **RISK 9** (9A), **RISK 5** (5A/5B) — sizing/overflow hardening.
6. **RISK 4** (4B), **RISK 7** (7A), **RISK 8** (8A), **RISK 11**,
   **ENHANCEMENTS E2–E4** — remaining hardening + cleanup.
7. **RISK 6** — last, pending the scoring-semantics decision.

Each numbered item = one (or a few small) commits, `cppcheck`-clean,
cross-compiled for ARM64/PPC64/Win64, with the swarm regression suite
green and the corresponding `pending_fixes.sh` test migrated into
`fixed_bugs.sh`.
