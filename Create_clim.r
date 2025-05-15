#!/bin/bash

# Define an array of month names
months=("Jan" "Feb" "March" "April" "May" "June" "July" "Aug" "Sep" "Oct" "Nov" "Dec")
module load CDO

# Define the base project directory
project_dir=""

# Loop over the months and execute the commands
for month in "${months[@]}"
do
    month_dir="${project_dir}/${month}"
    output_dir="${month_dir}/Complete"
    output_file="${output_dir}/pCO2_${month}.nc"

    # Navigate to the month's directory, check if it exists
    if [ -d "$month_dir" ]; then
        cd "$month_dir" || continue

        # Remove existing output file if it exists
        if [ -f "$output_file" ]; then
            rm -f "$output_file"
        fi

        # Compute the ensemble mean
        cdo ensmean *.nc "$output_file" && echo "Created ${output_file}"

        # Remove original NetCDF files after processing
        find "$month_dir" -maxdepth 1 -type f -name '*.nc' -delete && echo "Cleaned up ${month_dir}"

        # Return to the original directory
        cd - || exit
    else
        echo "Directory $month_dir does not exist; skipping $month."
    fi
done
