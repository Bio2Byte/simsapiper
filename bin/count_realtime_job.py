#This code sums the realtime of all jobs of the pipeline and gives it back in minutes
import csv
import sys

def extract_realtime_column(file_path):
    realtimes = []
    realtimes_tcoffee = []
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            realtimes.append(float(row['realtime']))
            for element in list(row.values()):
                if element.startswith("runTcoffee"):
                    realtimes_tcoffee.append(float(row['realtime']))
    return realtimes,realtimes_tcoffee

file_path = sys.argv[1]
realtime_columns,realtimes_tcoffee_columns = extract_realtime_column(file_path)
total = sum(realtime_columns)/60000 #convert to minutes
total_tcoffee = sum(realtimes_tcoffee_columns)/60000 #convert to minutes

print(f"Total realtime: {round(total)} minutes")
print(f"Total realtime tcoffee only: {round(total_tcoffee)} minutes")