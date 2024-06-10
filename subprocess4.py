import subprocess
import numpy as np
import random
import datetime
import psutil
import multiprocessing
import os
from concurrent.futures import ThreadPoolExecutor

cpu_cores = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

def generate_random_number():
    return random.randint(0, 10000)

# Define the compiled C++ program path
program_path = './dym-mod-metro'  # Use './program.exe' for Windows

# Beta1 values
#beta_0 = np.concatenate((np.linspace(0,1.8,120),np.linspace(1.85,3,round((3-1.85)/0.05))))
beta_0 = np.linspace(0,3,255)

# Run the compiled C++ program with parameters
def execute(parameters, ran=True):
    if ran == True:
        input = "./dym-mod-metro_n100 ./groups/myBT 4 4 5 0.000000e+00 " + str(float(parameters))+ " 0.000000e+00 " + str(generate_random_number())
        #input = ["./dym-mod-metro"] 
        #argument = ["./groups/myBO", "4", "4", "5", "0.000000e+00",str(float(parameters)),"0.000000e+00", str(generate_random_number())]
    if ran == False:
        input = "./dym-mod-metro ./groups/myBO 4 4 5 0.000000e+00 " + str(float(parameters))+ " 0.000000e+00 5137"
    #result = subprocess.run("taskset --cpu-list "+ str(core)+' '+input, shell=True, capture_output=True, text=True)
    process = subprocess.Popen(input, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
    #psutil.Process(process.pid).cpu_affinity([core])
    stdout, stderr = process.communicate()
    #return_code = process.wait()
    return f"Beta_1: {parameters} \n" + stdout
count = 0
print("Output:")
x = datetime.datetime.now()

value_tracker = 0
#num_processes = multiprocessing.cpu_count()
#pool = multiprocessing.Pool(processes=1)

#with open(f"./data/output{x}.txt", "a") as file:
num_instances = 3
def group_into_fives(arr):
    # Initialize an empty list to store the grouped subarrays
    grouped = []
    
    # Iterate over the array in steps of 5
    for i in range(0, len(arr), 15):
        # Append the subarray of the next 5 elements to the grouped list
        grouped.append(arr[i:i+15])
    
    return grouped
    
def run_cpp_instance(beta_val,core):
    os.sched_setaffinity(0, [core])
    process = subprocess.Popen("./dym-mod-metro_N300 ./groups/myBT 4 4 5 0.000000e+00 " + str(beta_val)+ " 0.000000e+00 " + str(generate_random_number()), stdout=subprocess.PIPE, shell=True)
    output, _ = process.communicate()
    result = output.decode("utf-8").strip()
    return beta_val, result
betas2 = group_into_fives(beta_0)
all_results = []
for set_v in betas2:
    count += 1
    print(f"Group calculation:  {count}")
    with ThreadPoolExecutor() as executor:
        # Submit tasks for each instance
        futures = [executor.submit(run_cpp_instance, i, core) for i, core in zip(set_v, cpu_cores)]

        # Gather results
        results = [future.result() for future in futures]
        executor.shutdown()
    all_results.append(results)
with open(f"./data/output{x}.txt", "a") as file:
    for res1 in all_results:
        for beta_val, result in res1:
            file.write(f"Beta_1: {beta_val} \n result: {result} \n")
#file.write(p1.result()+p2.result()+p3.result())
#pool.shutdown()
    




# Print the output and return code

