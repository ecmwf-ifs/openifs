#!/bin/bash

set -eu
set -o pipefail
#set -x

# Create Python 3 virtual environment for Ford
export FORDVENV="${TMPDIR}/ford"
python3 -m venv "${FORDVENV}"
source "${FORDVENV}/bin/activate"
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements-doc.txt

TEMP_PROJ_FILE="${PWD}/temp.md"
{
    cat ./ifs-test.md
    echo "revision: $(git rev-parse --abbrev-ref HEAD)"
    cat ../README.md
} > "${TEMP_PROJ_FILE}"

TEMP_FORT_FILE="${PWD}/temp.f90"
{
    echo "program empty"
    echo "! Temporary Fortran file to prevent Ford error when no Fortran files"
    echo "! are present"
    echo "end program empty"
} > "${TEMP_FORT_FILE}"

PAGES="${PWD}/pages"
mkdir -p "${PAGES}"
{
    echo "title: Test Documentation"
    echo
    echo "# Overview"
    echo
    echo "This section documents the test configurations."
} > "${PAGES}/index.md"

cd ../tests

# Keep the directory structure of ifs-tests/test
while IFS= read -r RESDIR ; do
    # Make the directory for the resolution
    mkdir -p "${PAGES}/${RESDIR}"
    # When making pages from Markdown Ford skips any directories that do
    # not contain index.md
    echo "title: $(basename "${RESDIR}") tests" > "${PAGES}/${RESDIR}/index.md"
done < <(find . -regextype posix-extended -type d -regex '\.\/t[lcoLCO]*[0-9]+')

# Make the main index.md for each test
while IFS= read -r TESTDIR ; do 
    mkdir -p "${PAGES}/${TESTDIR}"
    MDFILE="${PAGES}/${TESTDIR}/index.md" # Destination file
    {
        cat "${TESTDIR}/README.md"
        echo "# Parameters"
        echo
        echo "~~~~"
        {
            [[ -e "${TESTDIR}/params" ]] && cat "${TESTDIR}/params"
            grep -e NCONF -e TSTEP -e CUSTOP -e NPROMA "${TESTDIR}/setup" || true \
                | sed -E 's/[, ]*//g'
        } | sort
        echo "~~~~"
        echo
        echo "# Labels"
        echo
        echo "~~~~"
        grep "$(basename "${TESTDIR}")[ ]*APPEND PROPERTY" ./*/CMakeLists.txt \
            | tr -s " " \
            | cut -d" " -f6- \
            | sed 's/)//g' \
            | tr ' ' '\n' \
            | sort
        echo "~~~~"
    } >> "${MDFILE}"
done < <(find . -type d -name "test_*")

ford "${TEMP_PROJ_FILE}"

rm "${TEMP_PROJ_FILE}"
rm "${TEMP_FORT_FILE}"
rm -rf "${PAGES}"

exit 0
