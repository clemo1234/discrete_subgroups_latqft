import subprocess
import numpy as np
import random
import datetime
import psutil
import multiprocessing
from concurrent.futures import ThreadPoolExecutor

def generate_random_number():
    return random.randint(0, 10000)

# Define the compiled C++ program path
program_path = './dym-mod-metro'  # Use './program.exe' for Windows

# Beta1 values
#beta_0 = np.concatenate((np.linspace(0,1.8,120),np.linspace(1.85,3,round((3-1.85)/0.05))))
beta_0 = np.linspace(0,3,3)

# Run the compiled C++ program with parameters
def execute(parameters,ran=True):
    if ran == True:
        input = "./dym-mod-metro_n100 ./groups/myBO 4 4 5 0.000000e+00 " + str(float(parameters))+ " 0.000000e+00 " + str(generate_random_number())
        #input = ["./dym-mod-metro"] 
        #argument = ["./groups/myBO", "4", "4", "5", "0.000000e+00",str(float(parameters)),"0.000000e+00", str(generate_random_number())]
    if ran == False:
        input = "./dym-mod-metro ./groups/myBO 4 4 5 0.000000e+00 " + str(float(parameters))+ " 0.000000e+00 5137"
    result = subprocess.run(input, shell=True, capture_output=True, text=True)

    return "Beta_1: "+str(parameters)+"\n"+result.stdout
count = 0
print("Output:")
x = datetime.datetime.now()

value_tracker = 0
num_processes = multiprocessing.cpu_count()
#pool = multiprocessing.Pool(processes=1)
pool = ThreadPoolExecutor(3)
#with open(f"./data/output{x}.txt", "a") as file:
p1 = pool.submit(execute,beta_0[0])
p2 = pool.submit(execute,beta_0[1])
p3 = pool.submit(execute,beta_0[2])
print(p1.result())
#file.write(p1.result()+p2.result()+p3.result())
pool.shutdown()
    




# Print the output and return code

