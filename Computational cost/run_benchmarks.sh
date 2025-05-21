#!/bin/bash

# Ensure the output directory exists
mkdir -p "csv outputs"

# Compile all versions
echo "Compiling all versions..."

# Plain version
echo "Compiling plain versions..."
gfortran -cpp -DPRECISION_FP16 burgers_equation_supg.f90 -o burgers_supg_fp16
gfortran -cpp -DPRECISION_FP32 burgers_equation_supg.f90 -o burgers_supg_fp32
gfortran -cpp burgers_equation_supg.f90 -o burgers_supg_fp64
gfortran -cpp -DPRECISION_FP128 burgers_equation_supg.f90 -o burgers_supg_fp128

# Kahan summation version
echo "Compiling Kahan summation versions..."
gfortran -cpp -DPRECISION_FP16 burgers_equation_supg_ks.f90 -o burgers_supg_ks_fp16
gfortran -cpp -DPRECISION_FP32 burgers_equation_supg_ks.f90 -o burgers_supg_ks_fp32
gfortran -cpp burgers_equation_supg_ks.f90 -o burgers_supg_ks_fp64
gfortran -cpp -DPRECISION_FP128 burgers_equation_supg_ks.f90 -o burgers_supg_ks_fp128

# Create timing.csv if it doesn't exist
if [ ! -f timing.csv ]; then
    echo "precision,n_elements,version,seconds" > timing.csv
fi

# Loop over configurations and run tests
for n_elements in 32 64 128 256 512; do
    for precision in fp16 fp32 fp64 fp128; do
        for version in "plain" "ks"; do
            echo "Running test: n_elements=$n_elements, precision=$precision, version=$version"
            
            # Determine binary name
            if [ "$version" == "plain" ]; then
                binary="./burgers_supg_${precision}"
            else
                binary="./burgers_supg_ks_${precision}"
            fi
            
            # Run three times and get minimum time
            min_time=999999
            for i in 1 2 3; do
                echo "  Run $i/3..."
                
                # Use time command to measure execution time
                start_time=$(TIMEFORMAT='%3R'; { time $binary $n_elements > /dev/null 2>&1; } 2>&1)
                
                echo "    Time: $start_time seconds"
                
                # Convert to comparable format for bc
                start_time_bc=$(echo $start_time | tr -d 's')
                
                # Compare with min_time
                if (( $(echo "$start_time_bc < $min_time" | bc -l) )); then
                    min_time=$start_time_bc
                fi
            done
            
            # Append to timing.csv
            echo "$precision,$n_elements,$version,$min_time" >> timing.csv
            echo "  Best time: $min_time seconds"
        done
    done
done

echo "All tests completed. Results saved to timing.csv"

# 6666666666