#!/bin/bash

# Define the path to the Kraken2 database
DB_PATH="/home/people/s220868/project/databases/Bacteria_Kraken"

# Loop through each sample
# for i in {1..18}
# for i in $(seq 1 18)

start=1
end=18

for i in $(seq $start $end)
do
    echo "Processing Sample $i"
    
    # Define file names for paired and singleton reads
    R1="/home/people/s220868/project/sewage_data/reads/DTU2020-${i}-PRJ1066-RL-CPH-Sewage-${i}-af-5L_R1_001.trim.fq"
    R2="/home/people/s220868/project/sewage_data/reads/DTU2020-${i}-PRJ1066-RL-CPH-Sewage-${i}-af-5L_R2_001.trim.fq"
    singleton="/home/people/s220868/project/sewage_data/reads/DTU2020-${i}-PRJ1066-RL-CPH-Sewage-${i}-af-5L_R1_001.singletons.fq"

    # Create a directory for this sample's output
#    mkdir -p "./results/sample${i}"

    # Run kraken2 on paired reads
    kraken2 --threads 10 --db $DB_PATH --paired $R1 $R2 --output results_paired_${i}.txt --report report_paired_${i}.txt --report-zero-counts
    
    # Run kraken2 on singleton reads
    kraken2 --threads 10 --db $DB_PATH $singleton --output results_singletons_${i}.txt --report report_singletons_${i}.txt --report-zero-counts

    # Combine outputs for a comprehensive report
   # cat sample${i}/paired_output.txt sample${i}/singleton_output.txt > sample${i}/combined_output.txt
    
    # Generate report for the combined output
   # kraken2 --db $DB_PATH --report ./results/sample${i}/combined_report.txt --report-zero-counts < sample${i}/combined_output.txt
done

echo "All samples processed."
