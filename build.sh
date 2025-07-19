#!/bin/bash

# This script builds the IntNavLib library and all associated applications.
#
# Usage:
# ./build.sh [BuildType]
#
# BuildType can be "Debug" or "Release". It defaults to "Release".
#
# The script will:
# 1. Build IntNavLib with the specified build type.
# 2. Install it to $HOME/local_install.
# 3. Build all applications in the 'apps' directory, linking against the
#    just-installed IntNavLib.

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Configuration ---

# Determine build type from the first argument. Default to Release.
BUILD_TYPE="Release"
if [[ "$1" =~ ^([Dd]ebug)$ ]]; then
  BUILD_TYPE="Debug"
fi

# Set the installation prefix for the library.
# Using $HOME for a user-local installation.
INSTALL_PREFIX="$HOME/local_install"

# Get the absolute path of the script's directory (repository root).
REPO_ROOT=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

echo "================================================="
echo "IntNavLib Full Build Script"
echo "================================================="
echo "Build Type:         $BUILD_TYPE"
echo "Installation Prefix:  $INSTALL_PREFIX"
echo "Repository Root:      $REPO_ROOT"
echo "================================================="

# --- 1. Build and install IntNavLib library ---

echo
echo "--> Building and installing IntNavLib library..."
echo

# Create build directory for the library
LIB_BUILD_DIR="$REPO_ROOT/build"
mkdir -p "$LIB_BUILD_DIR"
cd "$LIB_BUILD_DIR"

# Configure with CMake
cmake .. -DCMAKE_BUILD_TYPE="$BUILD_TYPE" -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX"

# Build and install
make -j$(nproc)
make install

echo
echo "--> IntNavLib installed successfully to $INSTALL_PREFIX"
echo

# --- 2. Build all applications ---

echo "--> Building all applications in 'apps' directory..."
echo

APPS_DIR="$REPO_ROOT/apps"

for app_path in "$APPS_DIR"/*; do
  if [ -d "$app_path" ]; then
    APP_NAME=$(basename "$app_path")
    echo "-------------------------------------------------"
    echo "--> Building application: $APP_NAME"
    echo "-------------------------------------------------"

    if [ -f "$app_path/package.xml" ]; then
      echo "--> Detected ROS 2 package. Using colcon."
      cd "$app_path"
      colcon build --install-base "install" --cmake-args -DCMAKE_BUILD_TYPE="$BUILD_TYPE" -DCMAKE_PREFIX_PATH="$INSTALL_PREFIX"
    else
      echo "--> Detected standard CMake project."
      APP_BUILD_DIR="$app_path/build"
      mkdir -p "$APP_BUILD_DIR"
      cd "$APP_BUILD_DIR"
      cmake .. -DCMAKE_BUILD_TYPE="$BUILD_TYPE" -DCMAKE_PREFIX_PATH="$INSTALL_PREFIX"
      make -j$(nproc)
    fi
    echo "--> Finished building $APP_NAME."
    echo
  fi
done

echo "================================================="
echo "All builds completed successfully!"
echo "================================================="