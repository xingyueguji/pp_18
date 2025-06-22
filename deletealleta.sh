#!/bin/bash

# Store eta directories in a temporary file
tmpfile=$(mktemp)
find . -type d -name "eta" > "$tmpfile"

# Force input from terminal to ensure interactivity
while IFS= read -r dir; do
    echo "Found directory: $dir"
    echo -n "Do you want to delete '$dir' and all its contents? [y/N] "
    read confirm </dev/tty
    if [[ "$confirm" == [yY] ]]; then
        rm -rf "$dir"
        echo "Deleted $dir"
    else
        echo "Skipped $dir"
    fi
done < "$tmpfile"

rm "$tmpfile"