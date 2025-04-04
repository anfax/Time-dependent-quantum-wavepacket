import sys
import subprocess
subprocess.run(["make", "clean"])
subprocess.run(["make"])
# 检查是否有足够的命令行参数
if len(sys.argv) != 7:
    print("Usage: python script.py line_number position start_value end_value interval input_file")
    sys.exit(1)

# 获取命令行参数
line_number = int(sys.argv[1]) - 1  # 减1因为列表索引从0开始
position = int(sys.argv[2]) - 1  # 减1因为列表索引从0开始
start_value = int(sys.argv[3])
end_value = int(sys.argv[4])
interval = int(sys.argv[5])  # 新增的间隔参数
input_file = sys.argv[6]
import time 
# 遍历从start_value到end_value的数值，步长为interval
for value in range(start_value, end_value + 1, interval):
    # 读取文件内容
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # 分解原始行，得到数值列表
    original_line = lines[line_number]
    numbers = original_line.strip().split(',')

    # 替换第二个数值
    numbers[position] = str(value)

    # 组合新的数值列表为一行
    new_line = ','.join(numbers)

    # 写回新的行到文件
    with open(input_file, 'w') as file:
        lines[line_number] = new_line + '\n'
        file.writelines(lines)

    # 打印出替换后的行，以检查替换是否正确
    print(f"Line {line_number + 1} after replacement: {new_line}")

    # 执行run.py脚本
    subprocess.run(["python3", "run.py", str(value)])
    time.sleep(5.0)
