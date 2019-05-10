from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
import logging
log = logging.getLogger(__name__)
import numpy as np
import numpy.random as ran
import math, copy
import os
from itertools import repeat
from .distribution import get_value_from_distribution

np.set_printoptions(threshold=np.inf)

class Simulation(object):

    def __init__(self, configs):
        #read configs
        self.configs = configs

        #demographic parameters
        self.total_pop = configs["population"]
        self.labels = configs["labels"]

        self.infection_params = configs["infection"]
        self.infectivity = self.infection_params["infectivity"]

        #simulation parameters
        self.simulation = configs["simulation"]

        self.pop_labels = []
        self.pop_labels_by_name = {}
        self.pop_infection = np.zeros(self.total_pop)
        self.pop_infection_counter = np.zeros(self.total_pop)

        #list of population to exclude, e.g. recovered individuals
        self.pop_exclude = np.zeros(self.total_pop, dtype=bool)

        #use full matrix for now, but memory will blow out with too many individuals!
        self.pop_matrix = np.ones((self.total_pop, self.total_pop))

        #track new infections and infection network
        self.number_new_infections = 0
        self.infection_network = []

        #seed random number generator
        if "random_seed" in configs:
            ran.seed(configs["random_seed"])

    # @profile
    def populate(self):
        number = self.total_pop
        log.debug("populate: create {0} individuals.".format(number))

        #generate labels for every individual
        for i in range(number):
            lookup_table = {}
            log.debug("population id {0}:".format(i))
            for label_name in self.labels:
                if label_name not in self.pop_labels_by_name:
                    self.pop_labels_by_name[label_name] = []
                #generate label value
                label_value = get_value_from_distribution(self.labels[label_name])
                lookup_table[label_name] = label_value
                log.debug("pop_labels {0}: [{1}, {2}]".format(i, label_name, label_value))
                self.pop_labels_by_name[label_name].append(label_value)
            self.pop_labels.append(lookup_table)

        #generate population transmission matrix pop_matrix based on set up values
        for i in range(number):
            for j in range(number):
                if i == j:
                    self.pop_matrix[i,j] = 1
                else:
                    #base infectivity
                    self.pop_matrix[i,j] = self.infectivity

                    for label_name in self.labels:
                        transmission = self.labels[label_name]["transmission"]
                        if transmission["type"] == "bins":
                            bins_i = np.digitize([self.pop_labels[i][label_name]], transmission["bins"])[0] - 1
                            bins_j = np.digitize([self.pop_labels[j][label_name]], transmission["bins"])[0] - 1
                            self.pop_matrix[i, j] *= transmission["multiplier"][bins_i][bins_j]
                            #log.debug("individual {0}[{1}]: {2} in bins {3}".format(i,label_name, self.pop_labels[i][label_name], bins_i))
                            #log.debug("individual {0}[{1}]: {2} in bins {3}".format(j,label_name, self.pop_labels[j][label_name], bins_j))
                            #log.debug("matrix value is {0}".format(transmission["multiplier"][bins_i][bins_j]))
                        elif transmission["type"] == "list":
                            idx_i = transmission["list"].index(self.pop_labels[i][label_name])
                            idx_j = transmission["list"].index(self.pop_labels[j][label_name])
                            self.pop_matrix[i, j] *= transmission["multiplier"][idx_i][idx_j]
                        else:
                            print("error: unknown transmission type.")
                            exit(1)

        log.debug("pop_matrix is:\n{0}".format(self.pop_matrix))
        log.debug("pop_infection:{0}".format(self.pop_infection))


    def seed(self):
        num_infect = int(self.total_pop * self.infection_params["initial_prevalence"])
        log.info("seeding {0} infections.".format(num_infect))
        positions = np.random.randint(self.pop_infection.shape[0], size=num_infect)
        self.pop_infection[positions]=1
        for idx in range(num_infect):
            ran_duration = get_value_from_distribution(self.infection_params["duration_infection"])
            self.pop_infection_counter[positions[idx]] = ran_duration + 1
        log.debug("pop_infection_counter:\n{0}".format(self.pop_infection_counter))
        log.debug("pop_infection:{0}".format(self.pop_infection))


    # @profile
    def update(self,t, burn_in = False):
        #clear new infections and infection network
        self.number_new_infections = 0
        self.infection_network = []

        log.info("time = "+str(t))
        log.debug("pop_infection matrix at beginning:\n{0}".format(self.pop_infection))

        #step 1: find new infections based on interaction matrix
        #generate a random binomial matrix based on probabilities
        #log.debug("pop_infection:{0}".format(self.pop_infection))
        random_matrix_t = np.copy(self.pop_matrix)
        random_matrix_t[self.pop_infection > 0] = np.random.binomial(1,self.pop_matrix[self.pop_infection > 0])
        log.debug("transmission matrix after ran number:\n{0}".format(random_matrix_t))
        new_infection = np.dot(self.pop_infection, random_matrix_t)
        log.debug("new_infection matrix after dot:\n{0}".format(new_infection))
        new_infection[new_infection > 1]=1
        new_infection = new_infection * (1 - self.pop_exclude)
        new_infection_idx = ( new_infection - self.pop_infection ) > 0
        log.debug("new_infection_idx matrix:\n{0}".format(new_infection_idx))

        self.number_new_infections = np.sum(new_infection)

        #record infection network
        for i in range(self.total_pop):
            for j in range(self.total_pop):
                if i != j and self.pop_infection[j] != 1 and self.pop_exclude[j] != True:
                    if random_matrix_t[i][j] == 1:
                        #transmit from i to j
                        self.infection_network.append((i,j))

        for idx in range(len(new_infection_idx)):
            if new_infection_idx[idx] == True:
                ran_duration = get_value_from_distribution(self.infection_params["duration_infection"])
                self.pop_infection_counter[idx] = ran_duration + 1

        self.pop_infection = new_infection
        log.debug("pop_infection matrix after transmission:\n{0}".format(self.pop_infection))

        #step 2: check counter and clear infections
        active_infections = self.pop_infection_counter > 0
        new_recover_idx = self.pop_infection_counter == 1
        self.pop_infection_counter[active_infections] -= 1
        self.pop_exclude[new_recover_idx] = True

        log.debug("pop_infection_counter:\n{0}".format(self.pop_infection_counter))
        log.debug("pop_exclude:{0}".format(self.pop_exclude))
