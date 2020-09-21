import re
import os


def load(filename):
    # print("Loading data...")
    data = list(map(lambda x: x, open(str(filename)).readlines()))
    return data


def getcolumnames(data):
    columnames = re.sub('[#\[/\]]', '', data[0]).split()
    return columnames


def find_index(dummy_file, string):
    index_string = [string in x for x in dummy_file].index(True)
    return index_string


def make_and_cd(directory_name):
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
        print("Made " + directory_name)
    os.chdir(directory_name)
    print("Changed into directory: " + str(os.getcwd()))
