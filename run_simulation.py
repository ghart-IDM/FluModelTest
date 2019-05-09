#!/usr/bin/env python

from __future__ import print_function
from builtins import range
import yaml
import argparse, textwrap
import os
from datetime import datetime
import matplotlib.pyplot as plt

#### intro text ####

intro_text = textwrap.dedent('''\
=============================
Matrix-based Individual Model
Flu Prototype
=============================
''')

print(intro_text)

parser = argparse.ArgumentParser(
        prog='flu_prototype',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Flu Prototype")
parser.add_argument('--config', '-C', metavar='CONFIG_FILE_NAME', help="config file name", default="config.yaml")
parser.add_argument('--log', '-L', metavar='LOG_FILE_NAME', help="log file name", default="")
input_args = parser.parse_args()

#### setup logging ####
import logging

# #setup logging
log = logging.getLogger("flu")
log.setLevel(logging.INFO)

if input_args.log != "":
    ch = logging.FileHandler(input_args.log,mode='w')
else:
    ch = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s | %(name)-18s| %(levelname)-8s| %(message)s", "%H:%M:%S")
ch.setFormatter(formatter)
log.addHandler(ch)

from flu.simulation import Simulation
from flu.logs import setup_logging
from flu.resource import check_resources
from flu.report import Report

# @profile
def main():

    #read configs
    configs = yaml.load(open(input_args.config))

    #setup logging levels
    setup_logging("flu",configs["logging"])

    #check cpu and memory resources
    check_resources(configs)

    sim = Simulation(configs)
    report = Report(sim)

    #check output directory and create if it doesn't exist
    try:
        os.makedirs(configs["reports"]["output_directory"])
    except OSError:
        if not os.path.isdir(configs["reports"]["output_directory"]):
            raise

    if configs["reports"]["plots"]["enabled"] == True:
        #create the output directory
        #check plot directory and create if it doesn't exist
        try:
            os.makedirs(configs["reports"]["plots"]["output_directory"])
        except OSError:
            if not os.path.isdir(configs["reports"]["plots"]["output_directory"]):
                raise

    #simulation starts here

    log.info("simulation start.")

    sim.populate()

    seeding_time = configs["simulation"]["seeding_time"]

    for t in range(configs["simulation"]["total_time"]):

        if t in seeding_time:
            sim.seed()

        sim.update(t)
        report.update(t)

    log.info("simulation end.")

    #reports and plots
    report.write_reports()

    if configs["reports"]["plots"]["enabled"] == True:
        report.write_plots()

if __name__ == "__main__":
    main()
