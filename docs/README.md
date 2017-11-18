# monteising
montecarlo simulation of the ising model in root

Since the lattice is a mono dimensional `bool *` we need a method to
count the neighbours of each spin.

Here are examples of neighbour counting for mono dimensional and
bi dimensional lattice.

Reduntant terms were omitted: `(i // N) * N` is zero when in max
dimension, also, `(i + 1) % N` `(i - 1 -N) % N` may be reduntant.

We can easily seee a pattern. We use `i` as the index of the
spin in the lattice, `d` is a particular dimension where we
want to find the neighbours.
For a simple notation we used `//` and `**` in a python-ish style
meaning respectively integer division and power (power comes before
all the other operations).

```
i // N ** (d + 1)) * N ** (d + 1) + (i + N ** d         ) % N ** (d + 1)
i // N ** (d + 1)) * N ** (d + 1) + (i + N ** d - N ** d) % N ** (d + 1)
```
We can loop on `i` and `d` to obtain all the neighbours. The process
can be executed in parallel and be reduced with `if` statements on
the first and last iteration.

The first reduction is much more sensitive since let us to evaluate
a single power for each dimension > 1.

Multimimensional neighbour counting:

![1D Model][1Dmodel]

![2D model][2Dmodel]

[1Dmodel]: img/1D.png "1D Model"
[2Dmodel]: img/2D.png "2D Model"