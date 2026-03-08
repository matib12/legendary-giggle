#include <iostream>
#include <yaml-cpp/yaml.h>

int main() {
    YAML::Node config = YAML::LoadFile("simulation.yaml");

    double length = config["cavity"]["length"].as<double>();
    int finesse = config["cavity"]["finesse"].as<int>();
    double power = config["laser"]["power"].as<double>();

    std::cout << "Cavity length: " << length << std::endl;
    std::cout << "Laser power: " << power << std::endl;

    return 0;
}
