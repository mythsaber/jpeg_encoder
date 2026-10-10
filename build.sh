#!/usr/bin/env bash
set -e
CONFIG=${1:-Debug}
cd "$(dirname "$0")"
rm -rf build
cmake -S . -B build -A x64
cmake --build build --config "$CONFIG"
cmake --install build --config "$CONFIG"
