# Help

Use **'-h'** to show:

```
===== Center Protein Calc =====

Usage: [option] [value]
-h See this.
-i Specific input data path.
-r Specific reference data path.
-a Use 4 centrality algorithms together (BC, CC, DC, EC).
-b Use algorithm Betweenness Centrality (BC).
-c Use algorithm Closeness Centrality (CC).
-d Use algorithm Degree Centrality (DC).
-e Use algorithm Eigenvector Centrality (EC).

=========== Caution ===========

Must have at least one of ['-a', '-b', '-c', '-d', '-e'].
Must have '-r' and 'your refer data's path'.
Must have '-i' and  'your input data's path'.

============ Tips ============

Algorithm BC or CC will take a long trip to run (like O(N^3)), but it works!
Use '-b -c' together (save you 50% time)

============ About ============

Author: bipy@GitHub
Version: 20200707.4
```

**BC** and **CC** use "floyd" algorithm to calculate all the vertex, this precedure might take hours.

Luckily they share the same result of floyd, so use **'-b'** and **'-c'** together and gain a tea time : )

# Sample

Run DC and EC!

```shell
Center_Protein.exe -d -e -r "Refenence essential proteins.txt" -i "original dip.txt"
```

Run 4 algorithms together!

```shell
Center_Protein.exe -a -r "Refenence essential proteins.txt" -i "original dip.txt"
```

