#!/bin/bash
echo "Starting preprocessing…"
python post_upstream_processing.py metadata.csv data/quant/counts_data.txt data/quant/processed_count_data.csv
echo "Preprocessing completed."