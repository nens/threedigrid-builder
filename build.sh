#!/usr/bin/env bash
set -euo pipefail

BUILD_TYPE="Release"
BUILD_DIR="build"

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -d, --debug       Build in Debug mode (default: Release)"
    echo "  -r, --release     Build in Release mode (default)"
    echo "  -b, --build-dir   Set build directory (default: build)"
    echo "  -c, --clean       Clean build directory before building"
    echo "  -h, --help        Show this help message"
}

CLEAN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        -r|--release)
            BUILD_TYPE="Release"
            shift
            ;;
        -b|--build-dir)
            BUILD_DIR="$2"
            shift 2
            ;;
        -c|--clean)
            CLEAN=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Build type : ${BUILD_TYPE}"
echo "Build dir  : ${BUILD_DIR}"

if [[ $CLEAN -eq 1 && -d "${SCRIPT_DIR}/${BUILD_DIR}" ]]; then
    echo "Cleaning build directory..."
    rm -rf "${SCRIPT_DIR:?}/${BUILD_DIR}"
fi

cmake -S "${SCRIPT_DIR}" -B "${SCRIPT_DIR}/${BUILD_DIR}" -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
cmake --build "${SCRIPT_DIR}/${BUILD_DIR}" --config "${BUILD_TYPE}"

echo ""
echo "Copying extension module into threedigrid_builder/grid/fgrid/ ..."
cp "${SCRIPT_DIR}/${BUILD_DIR}"/_fgrid*.so "${SCRIPT_DIR}/threedigrid_builder/grid/fgrid/"

echo "Done."
