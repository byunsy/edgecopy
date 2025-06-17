import os
import sys
import time

BARLEN = 55

def make_title():
    title = """

     ____  ___   __    ____  __    ___   ___   _   
    | |_  | | \ / /`_ | |_  / /`  / / \ | |_) \ \_/
    |_|__ |_|_/ \_\_/ |_|__ \_\_, \_\_/ |_|    |_|  
    """
    try:
        _ = os.system("clear")
    except:
        _ = os.system("cls")
    print(title)


def make_header(text):
    time.sleep(0.25)
    print()
    print("="*BARLEN)
    print(text)
    print("="*BARLEN)
    print()


def print_command(argv):
    #print("python", end=" ")
    for arg in argv:
        if arg.startswith("-"):
            print(f"\n{arg}", end=" ")
        else:
            print(f"{arg}", end=" ")
    print()


def show_usage(prog_name):
    
    make_title()
    usage = '''
Usage:   exparascopy <command> <arguments>

Commands:

\033[32m[ Analyzing BAM files ]\033[0m
    depth       Calculates exonic read counts for all samples and for
                duplicated genes.

    agcn        Compute point and HMM aggregate copy number estimates.

    pscn        Compute allele-specific copy number estimates.

\033[32m[ General help ]\033[0m

    help        Show this help message.

'''
    print(usage)
