#!/bin/bash
# By Friedrich Rober
# go to root of repo
script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $script_dir/..
echo "Working in folder $(pwd)"

# Post-processing for the extracted examples from the documentation.
# - Move files into test_dir
test_dir="tst/Files/doc"
mkdir -p $test_dir
files=($(ls -1 tst/selfintersectingcomplexes*.tst))
echo "Found ${#files[@]} test file(s)"
for file in ${files[@]}; do
    echo "Processing $file"
    mv $file $test_dir/${file#"tst/"}
done