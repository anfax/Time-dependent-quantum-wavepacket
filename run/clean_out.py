import os 
#删除输出文件wf/* output/* 
def cleanoutput():
    for i in os.listdir('wf'):
        os.remove('wf/'+i)
    for i in os.listdir('output'):
        os.remove('output/'+i)
    
    
if __name__ == '__main__':
    cleanoutput()
    print('删除完成')
    #os.system('pause') 
    
