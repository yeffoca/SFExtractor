import data
import csv
import pandas as pd

for file in data.seCatsList:
    with open(file, 'r') as reader, open(f'{file[-15:-4]}CLEAN.csv', 'w') as writer:

        # Removing first 11 lines
        count = 0
        while count < 11:
            reader.readline()
            count += 1

        # writer = csv.writer(writer)
        # writer.writerow(data.header)
        rem_lines = reader.readlines()  # read all the remaining lines
        rem_lines_split = []
        for line in rem_lines:
            rem_lines_split.append(line.split())

        writer = csv.writer(writer)
        writer.writerow(data.header)
        writer.writerows(rem_lines_split)  # write all the remaining lines to another file
