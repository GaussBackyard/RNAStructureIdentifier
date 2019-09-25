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

## Prerequisites

To use the package, the only thing you need to install is Wolfram Mathematica [Free Trial](https://www.wolfram.com/mathematica/trial/).

We have tested our package using Mathematica versions 11.3.0 and 12.0.0 under UNIX and Windows system.

## Usage
Download the package to your Notebook directory and load the pacakge with
```
Get[NotebookDirectory[] <> "RNAStructureIdentifier.wl"]
```

## Quick Start

Usually you'll simply unpack the distribution tarball, configure and make:


## References

If you use our software package, you may want to cite the follwing publications:

[Thomas J.X. Li and Christian M. Reidys (2019)](http://arxiv.org/abs/1909.05744)
"On an enhancement of RNA probing data using Information Theory", arXiv:1909.05744


## Contact

We need your feedback! Send your comments, suggestions, and questions to
gauss.backyard@gmail.com

Thomas Li, Autumn 2019
