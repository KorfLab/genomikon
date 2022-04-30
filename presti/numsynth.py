import argparse
import numpy
import random

parser = argparse.ArgumentParser(
	description='Synthetic data generator for testing presti')
parser.add_argument('--count', type=int, metavar='<int>', required=True,
	help='number of regions to generate')
parser.add_argument('--depths', type=int, metavar='<ints>', required=True,
	nargs='+', help='list of classes/depths')
parser.add_argument('--widths', type=int, metavar='<ints>', required=True,
	nargs='+', help='list of peak widths')
parser.add_argument('--seed', type=int, metavar='<int>', required=False,
	help='use random seed')
arg = parser.parse_args()

pos = 1
for i in range(arg.count):
	depth = random.choice(arg.depths)
	width = random.choice(arg.widths)
	sd = depth ** 0.5
	print(f'# {pos} {pos+width-1} {depth}')
	for r in numpy.random.normal(depth, sd, width): print(f'{r:.1f}')
	pos += width

