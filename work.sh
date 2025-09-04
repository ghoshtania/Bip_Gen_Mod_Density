#!/bin/bash

# Input file is passed as the first argument
input_file=$1

# Count the number of edges (links)
edges=$(awk 'END{print NR}' "$input_file")
echo edges=$edges

# Get unique n1 nodes (red nodes)
awk '{print $1}' "$input_file" | sort -n | uniq > reds.txt
n1=$(awk 'END{print NR}' reds.txt)
echo reds=$n1

# Get unique n2 nodes (blue nodes)
awk '{print $2}' "$input_file" | sort -n | uniq > blues.txt
n2=$(awk 'END{print NR}' blues.txt)
echo blues=$n2

# Calculate degrees and weighted degrees for red nodes and save to a temporary file
awk '{deg[$1]++; weight[$1]+=$3} END {for (node in deg) print node, deg[node], weight[node]}' "$input_file" | sort -n > temp_red.txt

# Calculate degrees and weighted degrees for blue nodes and save to a temporary file
awk '{deg[$2]++; weight[$2]+=$3} END {for (node in deg) print node, deg[node], weight[node]}' "$input_file" | sort -n > temp_blue.txt

# Extract the second and third columns from temp_red.txt and save to degree_red.txt
awk '{print $2, $3}' temp_red.txt > degree_red.txt

# Extract the second and third columns from temp_blue.txt and save to degree_blue.txt
awk '{print $2, $3}' temp_blue.txt > degree_blue.txt
# Output the results to info.txt
echo "$n1 $n2 $edges" > info.txt

# Cleanup temporary files
rm blues.txt reds.txt temp_red.txt temp_blue.txt

# Create clean.txt with n1, n2, and weight from the original file
awk '{print $1, $2, $3}' "$input_file" > clean.txt
