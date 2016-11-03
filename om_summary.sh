#!/bin/sh

# Check how many lines there are in module - excludes blank lines
printf "\nNumber of lines of code & commets in OM: "
find ./om -name "*.py" -type f -exec grep . {} \; | wc -l
printf "\n"

# Check number of files using cloc
printf "\n CLOC OUTPUT: \n"
cloc om

# Run Tests
printf "\n RUN TESTS: \n"
py.test

# Find a way to get summary from pylint?
