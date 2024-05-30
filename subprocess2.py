import subprocess
import numpy as np
import random
import datetime
import psutil
import multiprocessing
import threading
from concurrent.futures import ThreadPoolExecutor

def generate_random_number():
    return random.randint(0, 10000)

# Define the compiled C++ program path
program_path = './dym-mod-metro'  # Use './program.exe' for Windows

# Beta1 values
#beta_0 = np.concatenate((np.linspace(0,1.8,120),np.linspace(1.85,3,round((3-1.85)/0.05))))
beta_0 = np.linspace(0,3,5)

# Run the compiled C++ program with parameters
def execute(parameters,core, ran=True):
    if ran == True:
        input = "./dym-mod-metro ./groups/myBO 4 4 5 0.000000e+00 " + str(float(parameters))+ " 0.000000e+00 " + str(generate_random_number())
        #input = ["./dym-mod-metro"] 
        #argument = ["./groups/myBO", "4", "4", "5", "0.000000e+00",str(float(parameters)),"0.000000e+00", str(generate_random_number())]
    if ran == False:
        input = "./dym-mod-metro ./groups/myBO 4 4 5 0.000000e+00 " + str(float(parameters))+ " 0.000000e+00 5137"
    #result = subprocess.run("taskset --cpu-list "+ str(core)+' '+input, shell=True, capture_output=True, text=True)
    process = subprocess.Popen(input, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
    psutil.Process(process.pid).cpu_affinity([core])
    stdout, stderr = process.communicate()
    return_code = process.wait()
    return f"Beta_1: {parameters} \n" + stdout
    #return result.stdout
count = 0
print("Output:")
x = datetime.datetime.now()
#with open(f"./data/output{x}.txt", "a") as file:
#    for value in beta_0:
#        count += 1
#        print(f"Step: {count}. Complete: {(count/len(beta_0))*100}%")
#        file.write(f"Beta_0: {value} \n" + execute(value,1))
value_tracker = 0
pool = ThreadPoolExecutor(5)
with open(f"./data/output{x}.txt", "a") as file:
    for value in range(5,len(beta_0)+1,5):
        count += 1
        betas_to_run = np.arange(value_tracker, value, step = 1)
        value_tracker += value
        print(f"Step: {count}. Complete: {(count/len(beta_0))*100}%")
        p1 = pool.submit(execute(beta_0[betas_to_run[0]],0))
        p2 = pool.submit(execute(beta_0[betas_to_run[1]],1))
        p3 = pool.submit(execute(beta_0[betas_to_run[2]],2))
        p4 = pool.submit(execute(beta_0[betas_to_run[3]],3))
        p5 = pool.submit(execute(beta_0[betas_to_run[4]],4))
        file.write(p1+p2+p3+p4+p5)    
        



# Print the output and return code

