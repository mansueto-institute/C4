## Algorithm Design

The aim of the software is to optimize compact, contiguous partitions of states.
Ultimately, efficient optimization 
  means efficiently calculating changes to compactness scores
  from the addition or removal of a single cell, rather than from scratch, each time.

For example, the IPQ (Polsby-Popper measure)
  relates the area of a shape to its perimeter.
It is trivial to add or remove one cell worth of area, 
  but keeping track of perimeters is a little harder.
Adding one cell,
  some parts of the old perimeter will now fall _within_ the shape,
  and those the same parts of the _new_ cell's own perimeter will too!
So one must track all of the adjacencies and add the right bits of cells'
  perimeters each time.
Every algorithm has considerations like this one.

The software is composed of three c++ classes.
These are exposed to python through cython.  

### Classes

#### Universe

This class is, in the districting context, the state.
The state contains a f

#### Region
#### Cells

### Contiguity Constraints

#### Connecting the Graph
#### Preserving Region Connectivity

### Optimization Routines

#### Greedy Optimizer
#### Split-Restart
#### Power 
#### Split-Line

### Objective Functions

### Cython Wrapping

### Python Script

Generating districts requires very specific geometries and topologies,
  and so these data are included with the software.
These data are loaded using python scripts that 
  users should NOT manipulate.

