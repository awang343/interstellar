#!/usr/bin/env bash
BIN="./projects_ray inputs/inis/antialias/required/sphere_checkerboard_bilinear.ini"

ls > test.log;

LOG_FILE="./combined.log"

echo "Logging to $LOG_FILE"
echo "Press Ctrl-C to stop."

while true; do
    timestamp=$(date +"%Y-%m-%d %H:%M:%S")

    echo "[$timestamp] === RUN START ==="
    $BIN > "$LOG_FILE" 2>&1
    exit_code=$?
    timestamp_end=$(date +"%Y-%m-%d %H:%M:%S")

    # If you want it to only re-run on nonzero exit, uncomment:
    if [[ ! $exit_code -eq 0 ]]; then break; fi

    sleep 1
done
