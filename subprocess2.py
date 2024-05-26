import subprocess
import numpy as np

# Define the compiled C++ program path
program_path = './dym-mod-metro'  # Use './program.exe' for Windows

# Beta0 values
beta_0 = np.linspace(0,5,50)


# Run the compiled C++ program with parameters
def execute(parameters):
    input = "./dym-mod-metro ./groups/myBT 4 4 4 " + str(float(parameters))+ " 0.000000e+00 0.000000e+00 5357"
    result = subprocess.run(input, shell=True, capture_output=True, text=True)
    return result.stdout
count = 0
print("Output:")
with open("output.txt", "a") as file:
    for value in beta_0:
        count += 1
        print(f"Step: {count}")
        file.write(f"Beta_0: {value} \n" + execute(value))
    



# Print the output and return code

