import argparse

parser = argparse.ArgumentParser(description='A test program.')

parser.add_argument("-p", "--print_string", help="Prints the supplied argument.", type=int)
parser.add_argument("-a", "--print_left", help="Prints the supplied left argument.", type=int)
args = parser.parse_args()

print(args.print_string, args.print_left)