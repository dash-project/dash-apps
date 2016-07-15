#!/bin/sh

BUILD_DIR=./build

FORCE_BUILD=false
if [ "$1" = "-f" ]; then
  FORCE_BUILD=true
fi

await_confirm() {
  if ! $FORCE_BUILD; then
    echo ""
    echo "   To build using these settings, hit ENTER"
    read confirm
  fi
}

mkdir -p $BUILD_DIR
rm -Rf $BUILD_DIR/*
(cd $BUILD_DIR && \
    cmake -DCMAKE_BUILD_TYPE=release \
          \
           ../ && \
    await_confirm && \
    make -j 4)
