# Help

Use **'-h'** to show:

```
===== Center Protein Calc =====

Usage: [option] [value]
-h See this.
-i Specific input data path.
-r Specific reference data path.
-s (optional) Specific step (default 100).
-a Use 5 centrality algorithms together (BC, CC, DC, EC, SigEP).
-b Use algorithm Betweenness Centrality (BC).
-c Use algorithm Closeness Centrality (CC).
-d Use algorithm Degree Centrality (DC).
-e Use algorithm Eigenvector Centrality (EC).
-p Use algorithm Significance-Based Essential Protein Discovery (SigEP).

=========== Caution ===========

Must have at least one of ['-a', '-b', '-c', '-d', '-e', '-p'].
Must have '-r' and 'your refer data's path'.
Must have '-i' and  'your input data's path'.

============ Tips ============

Algorithm BC or CC will take a long trip to run (like O(N^3)), but it works!
Use '-b -c' together (save you 50% time)

============ About ============

Author: bipy@GitHub
Version: 20200723.1
```

**BC** and **CC** use "floyd" algorithm to calculate all the vertex, this precedure might take some time.

Luckily they share the same result of floyd, so use **'-b'** and **'-c'** together and gain a tea time : )

# About SigEP

## Introduction

SigEP can identity the essential proteins from Protein-Protein Interaction network: we present a p-value calculation method for quantifying the statistical significance of each protein by considering both its degree and local clustering coefficient. To reduce the computational cost, we further present an upper bound of the p-value, which is less timeconsuming in practice. After calculating the p-value for each protein, we control the FDR of identified essential proteins using the well-known BH algorithm.

## Reference

Liu, Y., Liang, H., Zou, Q., & He, Z. (2020). Significance-Based Essential Protein Discovery. *IEEE/ACM Transactions on Computational Biology and Bioinformatics*.

## Calculate p_i

![](https://cdn.jsdelivr.net/gh/bipy/CDN@master/repo/Essential-Proteins/p.png)

![](https://cdn.jsdelivr.net/gh/bipy/CDN@master/repo/Essential-Proteins/c.png)

![](https://cdn.jsdelivr.net/gh/bipy/CDN@master/repo/Essential-Proteins/beta.png)



# Usage

Run DC and EC!

```shell
Center_Protein.exe -d -e -r "Reference essential proteins.txt" -i "original dip.txt"
```

Run 5 algorithms together!

```shell
Center_Protein.exe -a -r "Reference essential proteins.txt" -i "original dip.txt"
```

Run 5 algorithms together and specific step!

```shell
Center_Protein.exe -a -s 5 -r "Reference essential proteins.txt" -i "original dip.txt"
```

