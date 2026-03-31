# N-Body Gravitational Simulation

This project simulates gravitational interactions between particles using:

* Euler Method
* Runge-Kutta (RK4) Method

## Features

* Calculates gravitational acceleration using Newton's Law of Gravitation
* Tracks:

  * Momentum
  * Kinetic Energy
  * Gravitational Potential Energy
* Compares numerical methods (Euler vs RK4)
* Plots position, velocity, and acceleration over time

## Requirements

* Python 3
* numpy
* matplotlib

Install dependencies:

```
pip install numpy matplotlib
```

## How to Run

Make sure you have an `input.txt` file with particle data.

Then run:

```
python Simulation.py
```

## Output

* Energy and momentum comparisons printed in terminal
* Graphs displayed using matplotlib
* Data saved as `.npy` files

## Notes

* RK4 is more accurate than Euler for energy conservation
* Euler method may show drift over time

---

Created for physics simulation coursework.
