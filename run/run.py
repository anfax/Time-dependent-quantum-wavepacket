import os
import shutil
import subprocess
from datetime import datetime

# 获取脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))

# 获取上一级目录
parent_dir = os.path.dirname(script_dir)

# 提示用户输入新目录的名称
#new_dir_name = input("请输入新目录的名称: ")
# 从命令行参数获取改变的数值作为新目录的名称
import sys
if len(sys.argv) < 2:
    print("请提供数值作为新目录的名称")
    sys.exit(1)
new_dir_name = sys.argv[1]
# 获取当前时间并格式化
current_time = datetime.now().strftime("%Y%m%d%H%M%S")

# 创建新目录的完整路径，包括用户输入的名称和创建时间
new_dir_path = os.path.join(parent_dir, f"{new_dir_name}")#"_{current_time}")
os.makedirs(new_dir_path, exist_ok=True)

# 复制脚本所在目录下所有文件和文件夹到新创建的目录
for item in os.listdir(script_dir):
    src = os.path.join(script_dir, item)
    dst = os.path.join(new_dir_path, item)
    if os.path.isfile(src):
        shutil.copy2(src, dst)
    elif os.path.isdir(src):
        shutil.copytree(src, dst)

# 切换到新创建的目录
os.chdir(new_dir_path)
print('PWD:',new_dir_path)
import glob

if glob.glob('a*.out'):
    print('a*.out exist')
else:
    # 执行make命令
    subprocess.run(['make', 'clean'])
    make_result = subprocess.run(['make'], capture_output=True, text=True)
    print("Make output:")
    print(make_result.stdout)
    if make_result.returncode != 0:
        print("Make failed with error:")
        print(make_result.stderr)

# 执行submit命令
subprocess.run(['cp','sba.sh',current_time+'.sh'], capture_output=True, text=True)
submit_result = subprocess.run(['sbatch',current_time+'.sh'], capture_output=True, text=True)
print("Submit output:")
print(submit_result.stdout)
if submit_result.returncode != 0:
    print("Submit failed with error:")
    print(submit_result.stderr)
import time

# 查看squeue
squeue_result = subprocess.run(['squeue','-u','lih'], capture_output=True, text=True)
print("Squeue output:")
print(squeue_result.stdout)
if squeue_result.returncode != 0:
    print("Failed to get squeue output with error:")
    print(squeue_result.stderr)

