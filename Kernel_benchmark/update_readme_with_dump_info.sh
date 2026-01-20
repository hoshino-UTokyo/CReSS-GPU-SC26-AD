#!/bin/bash
# Update README.md files with dump data information

cd "$(dirname "$0")/.."

for readme in Kernel_benchmark/*/README.md; do
    kernel_dir=$(dirname "$readme")
    kernel_basename=$(basename "$kernel_dir")

    # Extract source file from README
    srcfile=$(grep -o "Src/[a-z0-9_]*\.f90" "$readme" | head -1)

    if [ -z "$srcfile" ] || [ ! -f "$srcfile" ]; then
        continue
    fi

    # Check if file has dump code
    if ! grep -q "call dump_init" "$srcfile"; then
        continue
    fi

    # Check if README already has Dump Data section
    if grep -q "## Dump Data Requirements" "$readme"; then
        continue
    fi

    # Extract dump information
    kernel_name=$(grep "call dump_init" "$srcfile" | head -1 | sed "s/.*dump_init('\([^']*\)').*/\1/")
    target=$(grep "DUMP_TARGET_" "$srcfile" | head -1 | sed 's/.*= *\([0-9]*\).*/\1/')

    # Create dump section
    dump_section=""
    dump_section+=$'\n## Dump Data Requirements\n\n'
    dump_section+="**Kernel Name**: \`$kernel_name\`\n"
    dump_section+="**Dump Target Call**: $target\n"
    dump_section+="**Dump Directory**: \`test_real/kernel_dump/$kernel_name/\`\n\n"

    dump_section+="### Input Scalars (params.txt)\n"
    dump_section+="\`\`\`\n"
    grep "call dump_scalar_" "$srcfile" | while read line; do
        varname=$(echo "$line" | sed "s/.*dump_scalar_[irc]('\([^']*\)'.*/\1/")
        echo "$varname"
    done | paste -sd ', ' >> /tmp/scalars_$$
    dump_section+="$(cat /tmp/scalars_$$ 2>/dev/null)\n"
    rm -f /tmp/scalars_$$
    dump_section+="\`\`\`\n\n"

    dump_section+="### Input Arrays\n"
    dump_section+="| File | Variable | Dimensions |\n"
    dump_section+="|------|----------|------------|\n"
    grep "call dump_array_" "$srcfile" | grep -v "_ref\." | while read line; do
        filename=$(echo "$line" | sed "s/.*('\([^']*\.bin\)'.*/\1/")
        varname=$(echo "$line" | sed "s/.*\.bin', *\([^,]*\),.*/\1/")
        dims=$(echo "$line" | sed "s/.*\.bin', *[^,]*, *\(.*\))/\1/" | tr -d ' ')
        echo "| \`$filename\` | \`$varname\` | \`$dims\` |"
    done >> /tmp/arrays_$$
    dump_section+="$(cat /tmp/arrays_$$ 2>/dev/null)\n\n"
    rm -f /tmp/arrays_$$

    dump_section+="### Output Arrays (Reference)\n"
    dump_section+="| File | Variable | Dimensions |\n"
    dump_section+="|------|----------|------------|\n"
    grep "call dump_array_" "$srcfile" | grep "_ref\." | while read line; do
        filename=$(echo "$line" | sed "s/.*('\([^']*\.bin\)'.*/\1/")
        varname=$(echo "$line" | sed "s/.*\.bin', *\([^,]*\),.*/\1/")
        dims=$(echo "$line" | sed "s/.*\.bin', *[^,]*, *\(.*\))/\1/" | tr -d ' ')
        echo "| \`$filename\` | \`$varname\` | \`$dims\` |"
    done >> /tmp/outarrays_$$
    dump_section+="$(cat /tmp/outarrays_$$ 2>/dev/null)\n"
    rm -f /tmp/outarrays_$$

    # Append to README
    echo "" >> "$readme"
    echo "$dump_section" >> "$readme"

    echo "Updated: $readme"
done
