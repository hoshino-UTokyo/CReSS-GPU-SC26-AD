#!/bin/bash
# Extract dump information from source files and generate dump_info.txt for each kernel

cd "$(dirname "$0")/.."

for srcfile in Src/*.f90; do
    basename=$(basename "$srcfile" .f90)
    
    # Check if file has dump code
    if ! grep -q "call dump_init" "$srcfile"; then
        continue
    fi
    
    # Extract kernel name from dump_init call
    kernel_name=$(grep "call dump_init" "$srcfile" | head -1 | sed "s/.*dump_init('\([^']*\)').*/\1/")
    
    if [ -z "$kernel_name" ]; then
        continue
    fi
    
    # Find matching kernel directory
    kernel_dir=$(ls -d Kernel_benchmark/*_${basename}* 2>/dev/null | head -1)
    if [ -z "$kernel_dir" ]; then
        kernel_dir=$(ls -d Kernel_benchmark/*_${kernel_name}* 2>/dev/null | head -1)
    fi
    
    if [ -z "$kernel_dir" ]; then
        continue
    fi
    
    # Extract dump calls and create dump_info.txt
    outfile="${kernel_dir}/dump_info.txt"
    
    echo "# Dump Data for kernel: $kernel_name" > "$outfile"
    echo "# Source: $srcfile" >> "$outfile"
    echo "" >> "$outfile"
    
    # Extract DUMP_TARGET value
    target=$(grep "DUMP_TARGET_" "$srcfile" | head -1 | sed 's/.*= *\([0-9]*\).*/\1/')
    echo "DUMP_TARGET: $target" >> "$outfile"
    echo "" >> "$outfile"
    
    echo "## Input Scalars (params.txt)" >> "$outfile"
    grep "call dump_scalar_" "$srcfile" | sed 's/^[ ]*/  /' >> "$outfile"
    echo "" >> "$outfile"
    
    echo "## Input Arrays" >> "$outfile"
    grep "call dump_array_" "$srcfile" | grep -v "_ref\." | sed 's/^[ ]*/  /' >> "$outfile"
    echo "" >> "$outfile"
    
    echo "## Output Arrays (Reference)" >> "$outfile"
    grep "call dump_array_" "$srcfile" | grep "_ref\." | sed 's/^[ ]*/  /' >> "$outfile"
    
    echo "Created: $outfile"
done
