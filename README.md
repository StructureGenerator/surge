<img src="resources/Ulogo.png" alt="drawing" width="100" align = "right"/>

# Surge: A Fast Open-Source Chemical Graph Generator
## About
Surge is a chemical structure generator based on the Nauty package and thereby on the principles of canonical augmentation.
More precisely, Surge generates all non-isomorphic constitutional isomers of a given molecular formula. [**Surge's article**](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-022-00604-9#citeas) is published in Journal of Cheminformatics and more details are given in its [**user manual**](https://github.com/StructureGenerator/SURGE/blob/main/doc/surge1_0.pdf) for more information. 

## Usage
Surge is a command line tool. Running `surge -u C10H16O` will generate the 452458 isomers of C<sub>10</sub>H<sub>16</sub>O in 0.1s on some vanilla flavor year-2021 PC. Running `surge -S C10H16O` outputs those structurs in SMILES format. You can either use `surge -S C10H16O > myresults.smi` to redirect the output into a result file, or use the `-o`switch to provide a filename. Further formats supported are SD Files (SDF) and a concise Surge-specific format.  
For large sets of structures, the -z option for compressing the output in gzip format will come in handy.

`surge -help` will show all options:

<p align="center">
  <img src="resources/logo.png" alt="drawing" width="400"
</p>

```
Usage: surge [-oFILE] [-z] [-u|-A|-S] [-T] [-e#|-e#:#] [-d#] [-c#] [-m#/#] formula

Make chemical graphs from a formula. Version 0.9.
Known elements are C,B,N,P,O,S,H,Cl,F,Br,I at their lowest valences.
Higher valences can be selected using Nx (Nitrogen/5), Sx,Sy (Sulfur 4/6), Px (Phosphorus/5).

formula = a formula like C8H6N2

  -E..  Define a new element (see the manual)
  -O#   Output stage: 1 after geng, 2 after vcolg, 3 after multig
        Default is to write SDF format
  -S    Output in SMILES format
  -A    Output in alphabetical format
  -u    Just count, don't write
  -e# -e#:#  Restrict to given range of distinct non-H bonds
  -t# -t#:#  Limit number of rings of length 3
  -f# -f#:#  Limit number of cycles of length 4
  -p# -p#:#  Limit number of cycles of length 5
  -b    Only rings of even length (same as only cycles of even length)
  -T    Disallow triple bonds
  -P    Require planarity
  -d#   Maximum degree not counting bond multiplicity or hydrogens (default 4)
  -c#   Maximum coordination number (default 4). This is the maximum number
        of distinct atoms (including H) that an atom can be bonded to
        Coordination number > 4 is only allowed if no neighbours are H
  -B#,...,# Specify sets of substructures to avoid (details in manual)
     1 = no triple bonds in rings up to length 7
     2 = Bredt's rule for two rings ij with one bond in
         common (33, 34, 35, 36, 44, 45)
     3 = Bredt's rule for two rings ij with two bonds in
         common (i,j up to 56)
     4 = Bredt's rule for two rings of length 6 sharing three bonds
     5 = no substructures A=A=A (in ring or not)
     6 = no substructures A=A=A in rings up to length 8
         For -B5 and -B6, the central atom only has 2 non-H neighbours
     7 = no K_33 or K_24 structure
     8 = none of cone of P4 or K4 with 3-ear
     9 = no atom in more than one ring of length 3 or 4
  -v     Write more information to stderr
  -m#/#  Do only a part. The two numbers are res/mod where 0<=res<mod.
  -oFILE Write the output to the given file rather than to stdout.
  -z     Write output in gzip format (only if compiled with zlib)
```
## Installation
### Option 1: Binary releases
Download one of the [releases](https://github.com/StructureGenerator/SURGE/releases) and run it.

### Option 2: Build from source code
1. Download the latest Nauty release from http://users.cecs.anu.edu.au/~bdm/nauty/ and build it following the instrucations on the page.
2. Download surge.c from the [source releases on this GitHub page](https://github.com/StructureGenerator/SURGE/releases) and put it into the nauty folder
3. Compile using the instructions at the beginning of surge.c. The following works on Linux, MacOS as well as with the https://MSYS2.org
Software Distribution and Building Platform for Windows. The latter was used to build the Windows release of Surge available on the [release page](https://github.com/StructureGenerator/SURGE/releases) .
```
gcc -o surge -O3 -DWORDSIZE=32 -DMAXN=WORDSIZE -DOUTPROC=surgeproc \
         -march=native -mtune=native -DPREPRUNE=surgepreprune \
         -DPRUNE=surgeprune -DGENG_MAIN=geng_main \
         surge.c geng.c planarity.c nautyW1.a
```
You can build-in gzip output using the zlib library (https://zlib.net). Add -DZLIB to the compilation, and link with the zlib library either by adding -lz or libz.a . This will activate the -z command to gzip the output :
```
gcc -o surge -O3 -DWORDSIZE=32 -DMAXN=WORDSIZE -DOUTPROC=surgeproc \
         -march=native -mtune=native -DPREPRUNE=surgepreprune -DZLIB \
         -DPRUNE=surgeprune -DGENG_MAIN=geng_main \
         surge.c geng.c planarity.c nautyW1.a -lz
```

### Option 3: Use docker 

First the container needs building:
```
docker build -t surge:latest .
```

Then it can be executed as follows:
```
docker run --rm -it surge:latest surge
```

## Misc
Surge was developed by [Brendan McKay](http://users.cecs.anu.edu.au/~bdm) with the help of [Christoph Steinbeck](https://github.com/steinbeck) and [Mehmet Aziz Yirik](https://github.com/mehmetazizyirik).
