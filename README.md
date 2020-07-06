# Help

Use **'-h'** to show:

```
==== Center Protein Calc ====
Usage: [option] [value]
-h See this.
-i Specific input data path.
-r Specific reference data path.
-d Use algorithm DC.
-c Use algorithm CC.
-b Use algorithm BC.
=========================
Author: bipy@GitHub
Version: 20200707.1
```

# Sample

```shell
Center_Protein.exe -d -r 'Refenence essential proteins.txt' -i 'original dip.txt'
```

# Caution

- Must have one of **'-b' '-c' '-d'**, the last one will be used if there is a dupe.
- Must have **'-r'** and **'your refer data's path'**.
- Must have **'-i'** and  **'your input data's path'**.
- Algorithm **BC** or **CC** will take a long trip to run (like $O(N^3)$), but it works!