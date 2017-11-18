# MonteIsing

## ToC
-[Overview](#Overview)
-[Lattice Class](#Lattice)
-[DrawLattice Class](#DrawLattice)
-[Multidimensional Neighbours Interactions and Counting](#MultiDim)

## Lattice Class

The class can build a multi dimensional square lattice of spins.

## DrawLattice Class

The class can draw 2D or 3D lattice from the class Lattice.

## Multimimensional Neighbour Interactions and Counting

Since the lattice is a mono dimensional `bool *` we need a method to
count the neighbours of each spin.

Here are examples of neighbour counting for mono dimensional and
bi dimensional lattice.

![1D Model][1Dmodel]

![2D model][2Dmodel]

Reduntant terms were omitted: `(i // N) * N` is zero when in max
dimension, also, `(i + 1) % N` `(i - 1 -N) % N` may be reduntant.

We can easily seee a pattern. We use `i` as the index of the
spin in the lattice, `d` is a particular dimension where we
want to find the neighbours.
For a simple notation we used `//` and `**` in a python-ish style
meaning respectively integer division and power (power comes before
all the other operations).

```
i // N**(d + 1)) * N**(d + 1) + (i + N**d             ) % N**(d + 1)
i // N**(d + 1)) * N**(d + 1) + (i - N**d + N**(d + 1)) % N**(d + 1)
```
We can loop on `i` and `d` to obtain all the neighbours. The process
can be executed in parallel and be reduced with `if` statements on
the first and last iteration.

The first reduction is much more sensitive since let us to evaluate
a single power for each dimension > 1.


[1Dmodel]: img/1D.png "1D Model"
[2Dmodel]: img/2D.png "2D Model"