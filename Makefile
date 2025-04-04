# 定义编译器
FC=ifort 

# 定义编译器 flags
FFLAGS= -O2 -fopenmp -qmkl

# 定义目标（最终的可执行文件）
TARGET=a.out

# 源文件目录
SRC_DIR=src app

# 对象文件目录
OBJ_DIR=obj

# 可执行文件目录
RUN_DIR=run

# 源文件列表，使用 wildcard 函数自动查找所有 .f 和 .f90 文件
SRCS=$(foreach dir,$(SRC_DIR),$(wildcard $(dir)/*.f) $(wildcard $(dir)/*.f90))

# 对象文件列表，将 .f 和 .f90 替换为 .o，并添加路径
OBJS=$(patsubst %.f,$(OBJ_DIR)/%.o,$(patsubst %.f90,$(OBJ_DIR)/%.o,$(notdir $(SRCS))))

# 默认目标是构建可执行文件
all: $(RUN_DIR)/$(TARGET)

# 创建必要的目录
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(RUN_DIR):
	mkdir -p $(RUN_DIR)

# 构建可执行文件
$(RUN_DIR)/$(TARGET): $(OBJS) | $(RUN_DIR)
	$(FC) $(FFLAGS) $^ -o $@

# 编译Fortran源文件，添加 -module 选项指定 .mod 输出目录
$(OBJ_DIR)/%.o: %.f | $(OBJ_DIR)
	$(FC) $(FFLAGS) -module $(OBJ_DIR) -c $< -o $@

$(OBJ_DIR)/%.o: %.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -module $(OBJ_DIR) -c $< -o $@

# 查找源文件的路径
vpath %.f $(SRC_DIR)
vpath %.f90 $(SRC_DIR)

# 清理目标，删除生成的文件（包括 .mod）
clean:
	rm -rf $(OBJ_DIR) 

# phony targets
.PHONY: all clean