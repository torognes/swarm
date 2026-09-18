#!/bin/sh
# Print swarm's version, after checking that every place carrying it agrees.
#
#   usage: version.sh [tag]
#
# swarm has no configure.ac to act as a single source of truth, so the
# version lives in four places that a release has to keep in step:
#
#   src/cli.cpp    the string "swarm -v" prints
#   man/swarm.1    the .TH header
#   CITATION.cff   what Zenodo and citation managers read
#   the git tag    what the release page is built from
#
# Run from the top of the source tree. The tag is only checked when given;
# both "v3.2.0" and "3.2.0" are accepted. The version is written to
# GITHUB_OUTPUT when that variable is set, so that the workflow can hand it
# to every other job.
set -eu

fail() { echo "::error::$*" >&2 ; exit 1 ; }

cli=$(sed -n 's/.*swarm_version *{ *"\([^"]*\)" *}.*/\1/p' src/cli.cpp)
man=$(sed -n 's/^\.TH .*"version \([^"]*\)".*/\1/p' man/swarm.1)
cff=$(sed -n 's/^version: *\([^ ]*\).*/\1/p' CITATION.cff | tr -d '\r')

test -n "${cli}" || fail "cannot read swarm_version from src/cli.cpp"
test -n "${man}" || fail "cannot read the version from the .TH header of man/swarm.1"
test -n "${cff}" || fail "cannot read the version from CITATION.cff"

echo "src/cli.cpp:  ${cli}" >&2
echo "man/swarm.1:  ${man}" >&2
echo "CITATION.cff: ${cff}" >&2

test "${cli}" = "${man}" || fail "man/swarm.1 says ${man}, src/cli.cpp says ${cli}"
test "${cli}" = "${cff}" || fail "CITATION.cff says ${cff}, src/cli.cpp says ${cli}"

# Accept the tag with or without its leading "v": releases are tagged
# v3.2.0, but a manual run may well be given the bare version.
tag="${1:-}"
if [ -n "${tag}" ] ; then
  echo "tag:          ${tag}" >&2
  test "${tag#v}" = "${cli}" || fail "tag ${tag} does not match version ${cli}"
fi

test -z "${GITHUB_OUTPUT:-}" || echo "version=${cli}" >> "${GITHUB_OUTPUT}"
echo "${cli}"
