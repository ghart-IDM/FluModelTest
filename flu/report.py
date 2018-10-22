import logging
log = logging.getLogger(__name__)
import numpy as np
import csv, json
import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

class Report:
    def __init__(self, sim):
        self.sim = sim
        self.configs = sim.configs
        self.new_infections = {}
        self.infection_label_counts = {}
        self.infection_networks = {}

    def update(self, t):
        self.report_infection(t)

    #reports and plots
    def write_reports(self):
        self.write_infection(self.configs["reports"]["output_directory"])
        self.write_infection_network(self.configs["reports"]["output_directory"] + "/infection.csv")

    def write_plots(self):
        self.plot_infection_network(self.configs["reports"]["plots"]["output_directory"] + "/infection_network.pdf")
#        self.plot_infection_label(self.configs["reports"]["plots"]["output_directory"] + "/infection_label.pdf")

    def report_infection(self,t):
        self.new_infections[t] = self.sim.number_new_infections
        self.infection_networks[t] = self.sim.infection_network

        #clear count by label
        count_by_label = {}
        for ind_label in self.sim.pop_labels:
            idx_indlabel = self.sim.pop_labels.index(ind_label)
            for label_name in ind_label:
                if label_name not in count_by_label:
                    count_by_label[label_name] = {}

                if self.configs["reports"]["labels"][label_name]["type"] == "bins":
                    idx = int(np.digitize(ind_label[label_name], self.configs["reports"]["labels"][label_name]["values"]))
                    bin_name = str(self.configs["reports"]["labels"][label_name]["values"][idx-1])+"-"+str(self.configs["reports"]["labels"][label_name]["values"][idx])
                elif self.configs["reports"]["labels"][label_name]["type"] == "list":
                    bin_name = ind_label[label_name]

                if bin_name not in count_by_label[label_name]:
                    count_by_label[label_name][bin_name] = 0

                if self.sim.pop_infection[idx_indlabel] == True:
                    count_by_label[label_name][bin_name] += 1

        for label_name in count_by_label:
            if label_name not in self.infection_label_counts:
                self.infection_label_counts[label_name] = {}
            for bin_name in count_by_label[label_name]:
                if bin_name not in self.infection_label_counts[label_name]:
                    self.infection_label_counts[label_name][bin_name] = []

                self.infection_label_counts[label_name][bin_name].append(count_by_label[label_name][bin_name])

    def write_infection(self, output_file_dir):
        for label_name in self.infection_label_counts:
            with open(output_file_dir+"infection_"+label_name+".csv", "wb") as outfile:
                writer = csv.writer(outfile)
                writer.writerow(self.infection_label_counts[label_name].keys())
                writer.writerows(zip(*self.infection_label_counts[label_name].values()))


    def write_infection_network(self, output_file_name):
        with open(output_file_name, 'wb') as outfile:
            w = csv.writer(outfile, delimiter=',')
            w.writerow(['time', 'from', 'to'])

            for sim_time,  networks in self.infection_networks.items():
                for from_id, to_id in networks:
                    w.writerow([sim_time, from_id, to_id])

    def plot_infection_network(self, output_file_name):
        G = nx.DiGraph()
        for sim_time,  networks in self.infection_networks.items():
            for from_id, to_id in networks:
                G.add_edge(from_id, to_id)

        pos = nx.layout.spring_layout(G)
        nodes = nx.draw_networkx_nodes(G, pos, node_color='yellow')
        edges = nx.draw_networkx_edges(G, pos, arrowstyle='->', arrowsize=10, width=2)
        nx.draw_networkx_labels(G, pos, font_size=10, font_family='sans-serif')
        plt.axis('off')
        plt.show()
