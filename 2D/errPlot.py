import pandas as pd
import matplotlib.pyplot as plt

# Step 1: Read the file
file_path = 'results/error.txt'
data = pd.read_csv(file_path)
data.columns = data.columns.str.strip()

# Step 2: Plot the data
plt.figure(figsize=(10, 6))
plt.plot(data['Iteration'], data['Error'], marker='o', linestyle='-', color='b')

# Step 3: Customize the plot
plt.title('Error vs. Iteration')
plt.xlabel('Iteration')
plt.ylabel('Error')
plt.grid(True)

# Step 4: Display the plot
plt.show()
