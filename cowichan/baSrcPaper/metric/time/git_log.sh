#!/bin/sh
git log --pretty=format:"%ad %s" > log.txt
tac log.txt > log_reverse.txt
rm log.txt
