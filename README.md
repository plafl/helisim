# Introduction

This repository contains a dump of my Master's thesis in Aerospace
Engineering that I finished in 2006, an helicopter simulator.

This is abandonware but I thought it might prove useful for someone. I
personally don't have any plans to update and maintain it, at least in
the short and mid term (unless paid :) ).

If you can read Spanish have a look at the [thesis itself](doc/memoria.pdf).

Code in `src` is under MIT license (see `src/LICENSE`).

Documentation in `doc` is under Creative Commons Attribution 4.0 (see `doc/LICENSE`).

# Features

The bad:

- It's written in Python 2, using old libraries.

- I didn't have much experience back then.

- The ground collision model sucks completely. Just throw it away. Who
  would have know after fighting with aerodynamics that rigid bodies
  still had surprises?

- Documentation is written in Spanish.

- Uses FlightGear for the graphical interface. I doubt it still works.

- It's just a flight dynamics model: no instrumentation, no fuel
  consumption, just the bare physics.

The good:

- Blade flapping.

- Finite state wake model plus blade element model for rotor
  aerodynamics.

- Wake effect on tail rotor.

- A simple ground effect model.

# Installation

Hahahaha!

