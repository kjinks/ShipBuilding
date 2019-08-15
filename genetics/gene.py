"""
Ken Jinks May 2018
file: gene.py
This file serves the base for functions and classes that deal with genetic
algorithms.
"""

# import numpy as np
from itertools import product
import math
import random
import io
import copy

def lerp(v, d):
    """
    linear interpolation
    """
    return v[0] * (1 - d) + v[1] * d

class GeneticSequence():
    """
    models a genetic sequence of n-genes
    that chromosomes can be read from as normalized values
    """
    chromosome_length = 4
    chromosome_index = 0
    index_stack = []

    def __init__(self, num_genes=4096, chromosome_length=4, filename=None):
        """
        initialize the sequence with a random set of bytes num_genes long
        """
        if filename == None:
            self.sequence = bytearray(random.randint(0,255) for _ in range(num_genes))
        else:
            self.load(filename)
        self.num_genes = len(self.sequence)
        self.chromosome_length = chromosome_length #in genes long

    def save(self, filename=str(random.randint(0,999999))):
        save_file = io.open(filename, 'wb')
        save_file.write(self.sequence)
        save_file.close()

    def load(self, filename=None):
        ##up to you to ensure correct chromosome length and gene size
        open_file = io.open(filename, 'rb')
        self.sequence = bytearray(open_file.read())
        open_file.close()
        self.reset()

    def reset(self):
        self.chromosome_index = 0
        self.index_stack = []
        self.index_stack.append(self.chromosome_index)

    def get_chromosome(self, index=0):
        """
        returns the array of genes that make up a chromosome
        """
        index = index % self.chromosome_length

        chromosome = (self.sequence[index:] + self.sequence[:index])[:self.chromosome_length]
        # print(self.chromosome_index, ' getting chromosome ', chromosome[0], chromosome[1])

        return chromosome

    def read_chromosome_value(self, chromosome, min=0.0, max=1.0):
        """
        converts the chromosome to a normalized
        value between min and max inclusively
        """
        total = 0.0


        min = float(min)
        max = float(max)

        for i in range(len(chromosome)):
            gene = self.sequence[(self.chromosome_index + i) % self.num_genes]
            total += gene * (256.0 ** i)

        total = total / (256.0 ** len(chromosome))

        #linear interpolation between min and max of total
        total = float(lerp([min, max], total))

        return total

    def read_next_value(self, min=0.0, max=1.0):
        """
        reads the chromosome that is currently pointed to
        and increments the chromosome index
        """
        value = self.read_chromosome_value(
            chromosome = self.get_chromosome(self.chromosome_index),
            min = min,
            max = max
        )
        # print(self.chromosome_index, ' read next value ', value, ' min ', min, ' max ', max)
        self.chromosome_index = (self.chromosome_index + self.chromosome_length) % self.num_genes
        # print('chromosome length ',self.chromosome_length)
        return value

    def index_push(self):
        self.index_stack.append(self.chromosome_index)
        # print('pushed_index {} stack size {}'.format(self.chromosome_index, len(self.index_stack)))

    def index_pop(self):
        self.chromosome_index = self.index_stack.pop()

        # print('popped_index {} stack size {}'.format(self.chromosome_index, len(self.index_stack)))

    def mutate(self, mutation_factor=0.1):

        for b in range(len(self.sequence)):
            if random.random() < mutation_factor:
                self.sequence[b] = (self.sequence[b] + random.randint(0, 255)) % 255

    def breed(self, donor):
        child = copy.copy(donor)
        for i in range(len(child.sequence)):
            if random.random() < 0.5:
                child.sequence[i] = self.sequence[i % len(self.sequence)]
        return child

    @staticmethod
    def create_test_data(gs):
        test_data = []
        for i in range(10):
            gs.index_push()
            test_data = test_data + [gs.read_next_value()]
            for j in range(10):
                gs.index_push()
                test_data = test_data + [gs.read_next_value()]
                for k in range(10):
                    gs.index_push()
                    test_data = test_data + [gs.read_next_value()]
                    gs.index_pop()
                gs.index_pop()
            gs.index_pop()
        return test_data

    @staticmethod
    def create_test_file(filename='test.tst'):
        gs = GeneticSequence()

        gs.save(filename="./test_000")

        result_file = io.open(filename, 'w')

        for i in GeneticSequence.create_test_data(gs):
            result_file.write("{}\n".format(i))

        result_file.close()

    @staticmethod
    def analyze_test_file(filename='./test.tst'):
        gs = GeneticSequence()

        gs.load(filename="./test_000")
        total_lines = 0

        test_data = GeneticSequence.create_test_data(gs)
        print('test_data has {} many points'.format(len(test_data)))
        result_file = io.open(filename, 'r')
        for line in result_file:
            expected = test_data.pop(0)
            result = float(line)
            if expected != result:
                print('expected {} result {}'.format(expected, result))
            else:
                total_lines += 1

        result_file.close()
        print("{} lines passed test".format(total_lines))
#gs = GeneticSequence()
#print (gs.get_chromosome(99))
#print (gs.read_chromosome_value(gs.get_chromosome(gene_length=8), min=-100.0, max=100.0))
