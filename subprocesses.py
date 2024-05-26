import subprocess

# Define the compiled C++ program path
program_path = './dym-mod-metro ./groups/myBT 3 4 4 0.000000e+00 0.000000e+00 0.000000e+00 5357'  # Change this to the path of your compiled C++ file

# Define the list of inputs
inputs = [
    "./groups/myBT 3 4 4 0.000000e+00 0.000000e+00 0.000000e+00 5357",
    "input2",
    "input3"
]

# Function to run the compiled C++ program with a specific input
def run_program(input_data):
    try:
        # Execute the compiled C++ program with the given input
        result = subprocess.run(program_path, input=input_data, text=True, capture_output=True, check=True)
        
        # Capture the output and return code
        output = result.stdout
        return_code = result.returncode
        
        # Print the output and return code
        print(f"Input: {input_data}")
        print(f"Output: {output}")
        print(f"Return Code: {return_code}")
        print('-' * 40)
        
    except subprocess.CalledProcessError as e:
        print(f"Error running program with input '{input_data}': {e}")

# Run the program with each input
for input_data in inputs:
    run_program(input_data)
