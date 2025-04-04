#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --partition all
##SBATCH -w xb001
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=Aoutput_%j.txt
#SBATCH --error=Aerror_%j.txt
#SBATCH --error=Aerror_%j.txt

export OMP_NUM_THREADS=16
source /etc/profile.d/modules.sh
module load compiler/2023.0.0
module load mkl
module load mpi
# module load conda
# conda activate
# module load gcc

# 记录开始时间
start_time=$(date +%s)
echo "Start Time: $(date -d @${start_time})"

# 执行您的程序
./a.out #python3 runtest.py 
# python3 name_of_program

# 记录结束时间
end_time=$(date +%s)
echo "End Time: $(date -d @${end_time})"

# 计算总用时
total_time=$(($end_time - $start_time))
echo "Total Time: $(($total_time / 60)) minutes and $(($total_time % 60)) seconds"

# 将开始时间、结束时间和总用时与路径一起记录
# echo "$(date -d @${start_time}) $(date -d @${end_time}) $(($total_time / 60)) minutes and $(($total_time % 60)) seconds $(pwd)" >>$HOME/program_end_paths_times.txt
echo "$(($total_time % 60)) seconds $(pwd)" >>$HOME/program_end_paths_times.txt

# 清理临时文件（可选）
rm -f start_time end_time total_time
