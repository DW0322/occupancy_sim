#!/usr/bin/env bash
# bootstrap.sh: install system deps, then configure, build and test.
#
# Linux and macOS only, and optional on both: it dispatches over
# apt/dnf/pacman/zypper/apk/brew and exits 1 when it finds none of them. The
# per-platform requirements, Windows included, are in README.md under
# Requirements; that list is the source of truth, this script is a shortcut.
#
# Needs cmake, a C++17 compiler, git, and the HDF5 DEVELOPMENT package (headers,
# not just the runtime); MorphIO won't configure without them.
#
# Installing system packages needs root. This script does NOT call sudo for you
# - set SUDO=sudo if your account needs it, which is the usual case:
#   SUDO=sudo ./bootstrap.sh
# Running as root already? Just ./bootstrap.sh.
#
# Once the deps are there you don't need this script again:
#   cmake --preset default && cmake --build --preset default && ctest --preset default
set -euo pipefail

if [ "${SUDO:-}" = "" ] && [ "$(id -u)" != "0" ]; then
  echo "note: not root and SUDO is unset - the package step will likely fail."
  echo "      re-run as:  SUDO=sudo $0"
  echo "      (harmless if the dependencies are already installed)"
fi

have() { command -v "$1" >/dev/null 2>&1; }

echo "== system packages"
if have apt-get;   then ${SUDO:-} apt-get update -qq && ${SUDO:-} apt-get install -y -qq cmake g++ git libhdf5-dev
elif have dnf;     then ${SUDO:-} dnf install -y cmake gcc-c++ git hdf5-devel
elif have pacman;  then ${SUDO:-} pacman -Sy --needed --noconfirm cmake gcc git hdf5
elif have zypper;  then ${SUDO:-} zypper --non-interactive install cmake gcc-c++ git hdf5-devel
elif have apk;     then ${SUDO:-} apk add --no-cache cmake g++ git hdf5-dev
elif have brew;    then brew install cmake git hdf5
else
  echo "!! unknown package manager. Install cmake, a C++17 compiler, git and"
  echo "!! the HDF5 development package, then re-run."
  exit 1
fi

echo "== python"
# installs the deps and the `occsim` command, so it works from any directory
pip install -q -e python 2>/dev/null \
  || pip install --break-system-packages -q -e python \
  || pip install --break-system-packages -q -r python/requirements.txt

echo "== configure + build"
cmake --preset default
cmake --build --preset default

echo "== tests"
ctest --preset default

cat <<'MSG'

built  build/occupancy_engine

If MorphIO was found already installed rather than fetched, you may need:
  export LD_LIBRARY_PATH=/usr/local/lib

smoke test:
  ./build/occupancy_engine --morph data/morph.asc --n 15 --div 15
  expect  = 81.7524 %   (pinned in tests/test_regression.cpp)

python side (occsim is a package - run it as a module, not as a file):
  cd python && python -m occsim run --n 15 --div 15
  or just `occsim run --n 15 --div 15` if the pip install above succeeded
MSG
