# RNAStructureIdentifier
A Mathematica package for identifying the target secondary structure from RNA structure ensemble.

Amongst other things, our implementations allow you to:

- compute the ensemble tree, a hierarchical bi-partition of the ensemble of structures
- calculate the base pairs having maximum information entropy
- predict the target structure via recursively querying about whether or not
	a base pair of maximum entropy is contained in the target
- visualize the ensemble tree, together with base pairs of maximum entropy
- draw a single secondary structure as a diagram
- present an ensemble of secondary structures as a greyscale diagram

## Introduction

The package focuses on how to identify a target from a Boltzmann ensemble of secondary structures. The key idea is to employ an information-theoretic approach to solve the problem, via considering a variant of the Rényi-Ulam game. Our framework is centered around the ensemble tree, a hierarchical bi-partition of the input ensemble, that is constructed by recursively querying about whether or not a base pair of maximum information entropy is contained in the target. These queries are answered via relating local with global chemical probing data, employing the modularity in RNA secondary structures (see References).

In [1](http://arxiv.org/abs/1909.05744), we present that  leaves  of the tree are comprised of sub-samples exhibiting a distinguished structure with high probability. In particular, for a Boltzmann ensemble incorporating probing data, which is well established in the literature, the probability of our framework correctly identifying the target in the leaf is greater than 90%.

## Installation

To use the package, the only thing you need to install is Wolfram Mathematica [Free Trial](https://www.wolfram.com/mathematica/trial/).

We have tested our package using Mathematica versions 11.3.0 and 12.0.0 under UNIX and Windows system.

## Usage
Download the package to your Notebook directory and load the pacakge with
```
Get[NotebookDirectory[] <> "RNAStructureIdentifier.wl"]
```

## Quick Start

We provide a quick tutorial (tutorial.nb) to demonstrate functions by examples.

## Functions

The Package includes the following functions:

| Program           | Description                                                                                                               |
| ------------------| :-------------------------------------------------------------------------------------------------------------------------|
| `Structure2Pair`  | Create a pair table from a dot-bracket notation of a secondary structure                            |
| `MaxEntfBP`       | Calculate the maximum entropy base pair for an ensemble of secondary structures.                    |
| `TreeMEBP`        | Calculate the ensemble tree, a hierarchical bi-partition for an ensemble of secondary structures    |
| `RanPdtMEBP`      | Calculate random walks on the ensemble tree based on the error rates                                |
| `DrawESB`         | Draw a greyscale diagram to present an ensemble of secondary structures                             |
| `DrawTree`        | Visualize the ensemble tree, together with base pairs of maximum entropy                            |

## References

If you use our package, you may want to cite the follwing publications:

1. [Thomas J.X. Li and Christian M. Reidys (2019)](http://arxiv.org/abs/1909.05744)
"On an enhancement of RNA probing data using Information Theory", arXiv:1909.05744


## Contact

We need your feedback! Send your comments, suggestions, and questions to
gauss.backyard@gmail.com

Thomas Li, Autumn 2019
