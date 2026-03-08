import yaml

with open("labs/YASML/python/simulation.yaml", "r") as f:
    config = yaml.safe_load(f)

length = config["cavity"]["length"]
finesse = config["cavity"]["finesse"]
power = config["laser"]["power"]

print("Cavity length:", length)
