import subprocess
import numpy as np
import random
import datetime

def generate_random_number():
    return random.randint(0, 100000)

# Define the compiled C++ program path
program_path = './dym-mod-metro'  # Use './program.exe' for Windows

# Beta0 values
beta_0 = np.linspace(0,1,50)


# Run the compiled C++ program with parameters
def execute(parameters, ran=True):
    if ran == True:
        input = "./dym-mod-metro ./groups/myBO 4 4 5 " + str(float(parameters))+ " 0.000000e+00 0.000000e+00 " + str(generate_random_number())
    if ran == False:
        input = "./dym-mod-metro ./groups/myBO 4 4 5 " + str(float(parameters))+ " 0.000000e+00 0.000000e+00 5137"
    result = subprocess.run(input, shell=True, capture_output=True, text=True)
    return result.stdout
count = 0
print("Output:")
x = datetime.datetime.now()
with open(f"./data/output{x}.txt", "a") as file:
    for value in beta_0:
        count += 1
        print(f"Step: {count}. Complete: {(count/len(beta_0))*100}%")
        file.write(f"Beta_0: {value} \n" + execute(value))
    



# Print the output and return code

