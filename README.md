# Help

Use **'-h'** to show:

```
===== Center Protein Calc =====

Usage: [option] [value]
-h See this.
-i Specific input data path.
-r Specific reference data path.
-d Use algorithm DC.
-c Use algorithm CC.
-b Use algorithm BC.

=========== Caution ===========

Must have at least one of '-b', '-c' or '-d'.
Must have '-r' and 'your refer data's path'.
Must have '-i' and  'your input data's path'.
Algorithm BC or CC will take a long trip to run (like O(N^3)), but it works!

============ Tips ============

Use '-b -c' together (save 50% time)

============ About ============

Author: bipy@GitHub
Version: 20200707.2

```

# Sample

```shell
Center_Protein.exe -b -c -r 'Refenence essential proteins.txt' -i 'original dip.txt'
```
