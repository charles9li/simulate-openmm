def run_simulation(simulation, simulation_options):
    if simulation_options.minimizeEnergy:
        simulation.minimizeEnergy()
    steps = simulation_options.steps
    simulation.step(steps)
