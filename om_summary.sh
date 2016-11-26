#!/bin/sh

# Check how many lines there are in module - excludes blank lines
printf "\n\n\n\n CHECK MODULE SIZE:"
printf "\nNumber of lines of code & comments in OM: "
find ./om -name "*.py" -type f -exec grep . {} \; | wc -l

# Check number of files using cloc
printf "\n\n\n CLOC OUTPUT (EXCLUDING TESTS): \n"
cloc om --exclude-dir='tests'

printf "\n\n\n CLOC OUTPUT - TEST FILES: \n"
cloc om/tests --exclude-dir='data'

# Prepare test data
printf "\n\n\n PREPARING TEST DATA: \n"
python om/tests/fake_data.py

# Run Tests & Check Coverage
printf "\n\n\n RUN TESTS & CHECK COVERAGE: \n"
coverage run --source om --omit="*/plts/*","*/fake_data.py" -m py.test
coverage report

# Find a way to get summary from pylint?

# Print out some new lines
printf "\n\n\n\n"
