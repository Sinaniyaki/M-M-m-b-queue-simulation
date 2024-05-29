# M/M/m/b Queue Simulation

## Overview
This project simulates an M/M/m/b queue, a queueing model where:
- **M**: Markovian (exponential) interarrival times
- **M**: Markovian (exponential) service times
- **m**: Multiple servers
- **b**: Finite buffer (queue size)

## Features
- Simulation of arrival and departure events
- Handling of multiple servers and finite queue
- Calculation of key metrics: average wait time, average queue length, server utilization, etc.
- Customizable parameters for arrival rate, service rate, number of servers, buffer size, number of departures, and random seed

## Getting Started

### Prerequisites
- C compiler (e.g., GCC)
- Standard C library

### Compilation
To compile the program, use the following command:

'''
gcc -o queue_simulation MMmb_queue.c -lm
'''

### Usage
Run the simulation with:

'''
./queue_simulation
'''

### Parameters
- `lambda`: Arrival rate
- `mu`: Service rate
- `m`: Number of servers
- `B`: Buffer size (maximum queue length)
- `D`: Number of departures to simulate
- `seed`: Random seed for reproducibility
