# Help

Use **'-h'** to show:

```
===== Center Protein Calc =====

Usage: [option] [value]
-h See this.
-i Specific input data path.
-r Specific reference data path.
-s (optional) Specific step (default 100).
-a Use 5 centrality algorithms together (BC, CC, DC, EC, DC_P).
-b Use algorithm Betweenness Centrality (BC).
-c Use algorithm Closeness Centrality (CC).
-d Use algorithm Degree Centrality (DC).
-e Use algorithm Eigenvector Centrality (EC).
-p Use algorithm Degree Centrality with p-value (DC_P).

=========== Caution ===========

Must have at least one of ['-a', '-b', '-c', '-d', '-e', '-p'].
Must have '-r' and 'your refer data's path'.
Must have '-i' and  'your input data's path'.

============ Tips ============

Algorithm BC or CC will take a long trip to run (like O(N^3)), but it works!
Use '-b -c' together (save you 50% time)

============ About ============

Author: bipy@GitHub
Version: 20200721.1
```

**BC** and **CC** use "floyd" algorithm to calculate all the vertex, this precedure might take some time.

Luckily they share the same result of floyd, so use **'-b'** and **'-c'** together and gain a tea time : )

# DC_P

![](https://cdn.jsdelivr.net/gh/bipy/CDN@master/repo/Essential-Proteins/dcp.png)

# Build

Require `gmp 6.20`



# Sample

Run DC and EC!

```shell
Center_Protein.exe -d -e -r "Refenence essential proteins.txt" -i "original dip.txt"
```

Run 5 algorithms together!

```shell
Center_Protein.exe -a -r "Refenence essential proteins.txt" -i "original dip.txt"
```

Run 5 algorithms together and specific step!

```shell
Center_Protein.exe -a -s 50 -r "Refenence essential proteins.txt" -i "original dip.txt"
```

