import argparse

from simulate.run import RunSimulation


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('parameter_file')
    command_line_args = parser.parse_args()

    simulation_information = RunSimulation(command_line_args.parameter_file)
    simulation_information.run()
