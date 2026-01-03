#!/usr/bin/env bash
set -euo pipefail

# ---- Step 1: rename *_A-B.txt -> *_A.0-B.0.txt (A,B are integers) ----
shopt -s nullglob

for f in *.txt; do
  # match "..._<int>-<int>.txt"
  if [[ "$f" =~ ^(.*_)([0-9]+)-([0-9]+)\.txt$ ]]; then
    pre="${BASH_REMATCH[1]}"
    a="${BASH_REMATCH[2]}"
    b="${BASH_REMATCH[3]}"
    new="${pre}${a}.0-${b}.0.txt"
    if [[ ! -e "$new" ]]; then
      mv -- "$f" "$new"
    fi
  fi
done

# ---- Step 2: copy requested bins ----
# helper: copy "prefix<old>.txt" -> "prefix<new>.txt" for all matching prefixes
copy_bin () {
  local old="$1"
  local new="$2"
  for src in *_"$old".txt; do
    [[ -e "$src" ]] || continue
    local dst="${src%_$old.txt}_$new.txt"
    if [[ ! -e "$dst" ]]; then
      cp -- "$src" "$dst"
    fi
  done
}

copy_bin "0.0-10.0"   "0.0-6.5"
copy_bin "10.0-20.0"  "6.5-14.5"
copy_bin "20.0-30.0"  "14.5-26.5"
copy_bin "30.0-100.0" "26.5-100.0"
copy_bin "0.0-100.0"  "0.0-15.0"
copy_bin "0.0-100.0"  "15.0-100.0"

echo "Done."
echo "Check counts per family, e.g.:"
echo "  ls PbPb_FA_tnpU_*.txt | wc -l"