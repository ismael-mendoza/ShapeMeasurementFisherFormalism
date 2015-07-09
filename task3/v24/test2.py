import csv


# import csv
# file1 = open('names.csv', 'r')
# reader = csv.reader(file1)
# for row in reader:
#     print row

with open('names.csv') as csvfile:
     reader = csv.DictReader(csvfile)
     for row in reader:
	     print(row['first_name'], row['last_name'])


# file2 = open('names.csv', 'w') ##append
# fieldnames = ['first_name', 'last_name']
# writer = csv.DictWriter(file2, fieldnames=fieldnames)
# writer.writeheader()
# writer.writerow({'first_name': 'Baked', 'last_name': 'Beans'})
# writer.writerow({'first_name': 'Lovely', 'last_name': 'Spam'})
# writer.writerow({'first_name': 'Wonderful', 'last_name': 'Spam'})