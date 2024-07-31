#!/bin/bash


# Clear previous output and prepare environment
rm -rf sewage.map
mkdir sewage.map


echo "RUNNING KMA ON SAMPLE 1"
kma -ipe reads/DTU2020-1-PRJ1066-RL-CPH-Sewage-1-af-5L_R1_001.trim.fq reads/DTU2020-1-PRJ1066-RL-CPH-Sewage-1-af-5L_R2_001.trim.fq -i reads/DTU2020-1-PRJ1066-RL-CPH-Sewage-1-af-5L_R1_001.singletons.fq -o sewage.map/sample1 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 2"
kma -ipe reads/DTU2020-2-PRJ1066-RL-CPH-Sewage-2-af-5L_R1_001.trim.fq reads/DTU2020-2-PRJ1066-RL-CPH-Sewage-2-af-5L_R2_001.trim.fq -i reads/DTU2020-2-PRJ1066-RL-CPH-Sewage-2-af-5L_R1_001.singletons.fq -o sewage.map/sample2 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 3"
kma -ipe reads/DTU2020-3-PRJ1066-RL-CPH-Sewage-3-af-5L_R1_001.trim.fq reads/DTU2020-3-PRJ1066-RL-CPH-Sewage-3-af-5L_R2_001.trim.fq -i reads/DTU2020-3-PRJ1066-RL-CPH-Sewage-3-af-5L_R1_001.singletons.fq -o sewage.map/sample3 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 4"
kma -ipe reads/DTU2020-4-PRJ1066-RL-CPH-Sewage-4-af-5L_R1_001.trim.fq reads/DTU2020-4-PRJ1066-RL-CPH-Sewage-4-af-5L_R2_001.trim.fq -i reads/DTU2020-4-PRJ1066-RL-CPH-Sewage-4-af-5L_R1_001.singletons.fq -o sewage.map/sample4 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 5"
kma -ipe reads/DTU2020-5-PRJ1066-RL-CPH-Sewage-5-af-5L_R1_001.trim.fq reads/DTU2020-5-PRJ1066-RL-CPH-Sewage-5-af-5L_R2_001.trim.fq -i reads/DTU2020-5-PRJ1066-RL-CPH-Sewage-5-af-5L_R1_001.singletons.fq -o sewage.map/sample5 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 6"
kma -ipe reads/DTU2020-6-PRJ1066-RL-CPH-Sewage-6-af-5L_R1_001.trim.fq reads/DTU2020-6-PRJ1066-RL-CPH-Sewage-6-af-5L_R2_001.trim.fq -i reads/DTU2020-6-PRJ1066-RL-CPH-Sewage-6-af-5L_R1_001.singletons.fq -o sewage.map/sample6 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 7"
kma -ipe reads/DTU2020-7-PRJ1066-RL-CPH-Sewage-7-af-5L_R1_001.trim.fq reads/DTU2020-7-PRJ1066-RL-CPH-Sewage-7-af-5L_R2_001.trim.fq -i reads/DTU2020-7-PRJ1066-RL-CPH-Sewage-7-af-5L_R1_001.singletons.fq -o sewage.map/sample7 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 8"
kma -ipe reads/DTU2020-8-PRJ1066-RL-CPH-Sewage-8-af-5L_R1_001.trim.fq reads/DTU2020-8-PRJ1066-RL-CPH-Sewage-8-af-5L_R2_001.trim.fq -i reads/DTU2020-8-PRJ1066-RL-CPH-Sewage-8-af-5L_R1_001.singletons.fq -o sewage.map/sample8 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 9"
kma -ipe reads/DTU2020-9-PRJ1066-RL-CPH-Sewage-9-af-5L_R1_001.trim.fq reads/DTU2020-9-PRJ1066-RL-CPH-Sewage-9-af-5L_R2_001.trim.fq -i reads/DTU2020-9-PRJ1066-RL-CPH-Sewage-9-af-5L_R1_001.singletons.fq -o sewage.map/sample9 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 10"
kma -ipe reads/DTU2020-10-PRJ1066-RL-CPH-Sewage-10-af-5L_R1_001.trim.fq reads/DTU2020-10-PRJ1066-RL-CPH-Sewage-10-af-5L_R2_001.trim.fq -i reads/DTU2020-10-PRJ1066-RL-CPH-Sewage-10-af-5L_R1_001.singletons.fq -o sewage.map/sample10 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 11"
kma -ipe reads/DTU2020-11-PRJ1066-RL-CPH-Sewage-11-af-5L_R1_001.trim.fq reads/DTU2020-11-PRJ1066-RL-CPH-Sewage-11-af-5L_R2_001.trim.fq -i reads/DTU2020-11-PRJ1066-RL-CPH-Sewage-11-af-5L_R1_001.singletons.fq -o sewage.map/sample11 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 12"
kma -ipe reads/DTU2020-12-PRJ1066-RL-CPH-Sewage-12-af-5L_R1_001.trim.fq reads/DTU2020-12-PRJ1066-RL-CPH-Sewage-12-af-5L_R2_001.trim.fq -i reads/DTU2020-12-PRJ1066-RL-CPH-Sewage-12-af-5L_R1_001.singletons.fq -o sewage.map/sample12 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 13"
kma -ipe reads/DTU2020-13-PRJ1066-RL-CPH-Sewage-13-af-5L_R1_001.trim.fq reads/DTU2020-13-PRJ1066-RL-CPH-Sewage-13-af-5L_R2_001.trim.fq -i reads/DTU2020-13-PRJ1066-RL-CPH-Sewage-13-af-5L_R1_001.singletons.fq -o sewage.map/sample13 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 14"
kma -ipe reads/DTU2020-14-PRJ1066-RL-CPH-Sewage-14-af-5L_R1_001.trim.fq reads/DTU2020-14-PRJ1066-RL-CPH-Sewage-14-af-5L_R2_001.trim.fq -i reads/DTU2020-14-PRJ1066-RL-CPH-Sewage-14-af-5L_R1_001.singletons.fq -o sewage.map/sample14 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 15"
kma -ipe reads/DTU2020-15-PRJ1066-RL-CPH-Sewage-15-af-5L_R1_001.trim.fq reads/DTU2020-15-PRJ1066-RL-CPH-Sewage-15-af-5L_R2_001.trim.fq -i reads/DTU2020-15-PRJ1066-RL-CPH-Sewage-15-af-5L_R1_001.singletons.fq -o sewage.map/sample15 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 16"
kma -ipe reads/DTU2020-16-PRJ1066-RL-CPH-Sewage-16-af-5L_R1_001.trim.fq reads/DTU2020-16-PRJ1066-RL-CPH-Sewage-16-af-5L_R2_001.trim.fq -i reads/DTU2020-16-PRJ1066-RL-CPH-Sewage-16-af-5L_R1_001.singletons.fq -o sewage.map/sample16 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 17"
kma -ipe reads/DTU2020-17-PRJ1066-RL-CPH-Sewage-17-af-5L_R1_001.trim.fq reads/DTU2020-17-PRJ1066-RL-CPH-Sewage-17-af-5L_R2_001.trim.fq -i reads/DTU2020-17-PRJ1066-RL-CPH-Sewage-17-af-5L_R1_001.singletons.fq -o sewage.map/sample17 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "RUNNING KMA ON SAMPLE 18"
kma -ipe reads/DTU2020-18-PRJ1066-RL-CPH-Sewage-18-af-5L_R1_001.trim.fq reads/DTU2020-18-PRJ1066-RL-CPH-Sewage-18-af-5L_R2_001.trim.fq -i reads/DTU2020-18-PRJ1066-RL-CPH-Sewage-18-af-5L_R1_001.singletons.fq -o sewage.map/sample18 -t_db /home/people/s220868/project/databases/bacdb/bacteria.ATG -ts 0 -oa -ef -nf -nc -tmp


echo "All files processed. KMA run complete."
