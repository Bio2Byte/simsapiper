#This code sums the realtime of all jobs of the pipeline and gives it back in minutes
import csv
import sys

def extract_realtime_column(file_path):
    realtimes = []
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            realtimes.append(float(row['realtime']))
    return realtimes

file_path = sys.argv[1]
realtime_columns = extract_realtime_column(file_path)
total = sum(realtime_columns)/60000 #convert to minutes
print(f"Total realtime: {round(total)} minutes")