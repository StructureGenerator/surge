<img src="resources/Ulogo.png" alt="drawing" width="100" align = "right"/>

# Surge: A Fast Open-Source Chemical Graph Generator
## About
Surge is a chemical structure generator based on the Nauty package and thereby on the principles of canonical augmentation.
More precisely, Surge generates all non-isomorphic constitutional isomers of a given molecular formula. [**Surge's article**](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-022-00604-9) is published in Journal of Cheminformatics and more details are given in its [**user manual**](https://github.com/StructureGenerator/SURGE/blob/main/doc/surge2_0.pdf).

## What's New in 2.0
- **Aromaticity filtering** (`-R`): Removes duplicate Kekule structures under carbon-ring aromaticity.
- **6-member cycle limits** (`-h#`, `-h#:#`): Restrict the number of cycles of length 6.
- **Chord-free 6-carbon cycle limits** (`-C#`, `-C#:#`): Restrict the number of chord-free cycles of 6 carbon atoms.
- **SDfile output flag** (`-F`): Explicit flag for SDfile output format.
- **64-bit word size**: Surge 2.0 is built with `WORDSIZE=64`, supporting larger molecules.
- Updated to nauty 2.9.3.

## Usage
Surge is a command line tool. Running `surge -u C10H16O` will generate the 452458 isomers of C<sub>10</sub>H<sub>16</sub>O in 0.1s on some vanilla flavor year-2021 PC. Running `surge -S C10H16O` outputs those structures in SMILES format. You can either use `surge -S C10H16O > myresults.smi` to redirect the output into a result file, or use the `-o` switch to provide a filename. Further formats supported are SD Files (SDF) and a concise Surge-specific format.
For large sets of structures, the -z option for compressing the output in gzip format will come in handy.

`surge -help` will show all options:

```
Usage: ./surge [-oFILE] [-z] [-A|-S|-F] [-T] [-e#|-e#:#] [-R] [-d#] [-c#] [-m#/#] formula

Make chemical graphs from a formula. Version 2.0.
  Known elements are C,B,N,P,O,S,H,Cl,F,Br,I at their lowest valences.
  Higher valences can be selected using Nx (Nitrogen/5), Sx,Sy (Sulfur 4/6)
   Px (Phosphorus/5).

  formula = a formula like C8H6N2

  -u    Just count, don't write molecules (default)
  -S    Output in SMILES format
  -F    Output in SDfile format
  -A    Output in alphabetical format
  -e# -e#:#  Limit the number of distinct non-H bonds
  -t# -t#:#  Limit the number of cycles of length 3
  -f# -f#:#  Limit the number of cycles of length 4
  -p# -p#:#  Limit the number of cycles of length 5
  -h# -h#:#  Limit the number of cycles of length 6
  -C# -C#:#  Limit the number of chord-free cycles 6 carbon atoms
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
  -R    Enable aromaticity detection (filters duplicate Kekule structures)
  -v    Write more information to stderr
  -m#/# Do only a part. The two numbers are res/mod where 0<=res<mod.
  -oFILE Write the output to the given file rather than to stdout.
  -E..  Define a new element (see the manual)
  -z    Write output in gzip format (only if compiled with zlib)
```
## Installation
### Option 1: Binary releases
Download one of the [releases](https://github.com/StructureGenerator/SURGE/releases) and run it. Pre-built binaries are available for Linux (x86_64), macOS (Apple Silicon), and Windows (x86_64). On Intel Macs, the Apple Silicon binary runs transparently via Rosetta 2, or you can build from source (see below).

### Option 2: Build from source code
1. Download the latest Nauty release from https://users.cecs.anu.edu.au/~bdm/nauty/ and build it following the instructions on the page.
2. Download the Surge source from the [source releases on this GitHub page](https://github.com/StructureGenerator/SURGE/releases) and put the `src/` directory alongside the nauty directory.
3. Compile using the Makefile in the `src/` directory. The following works on Linux, macOS, as well as with https://MSYS2.org on Windows:
```
cd src
make surge NAUTY=/path/to/nauty NAUTYLIB=/path/to/nauty
```
Or compile directly with gcc:
```
gcc -o surge -O3 -I /path/to/nauty -DWORDSIZE=64 -DMAXN=WORDSIZE \
         -DOUTPROC=surgeproc -march=native -DPREPRUNE=surgepreprune \
         -DPRUNE=surgeprune -DGENG_MAIN=geng_main \
         surge.c geng.c planarity.c /path/to/nautyL1.a
```
You can build-in gzip output using the zlib library (https://zlib.net). Add `-DZLIB` to the compilation, and link with the zlib library either by adding `-lz` or `libz.a`. This will activate the `-z` command to gzip the output:
```
gcc -o surge -O3 -I /path/to/nauty -DWORDSIZE=64 -DMAXN=WORDSIZE \
         -DOUTPROC=surgeproc -march=native -DPREPRUNE=surgepreprune \
         -DZLIB -DPRUNE=surgeprune -DGENG_MAIN=geng_main \
         surge.c geng.c planarity.c /path/to/nautyL1.a -lz
```

**Note for macOS Homebrew users:** If you installed nauty via Homebrew (`brew install nauty`), the library file may be named `libnautyL1.a` rather than `nautyL1.a`. You can create a symlink to make it work with the Makefile:
```
NAUTY_LIB=$(brew --prefix nauty)/lib
ln -s "$NAUTY_LIB/libnautyL1.a" "$NAUTY_LIB/nautyL1.a"
make surge NAUTY=$(brew --prefix nauty)/include/nauty NAUTYLIB=$(brew --prefix nauty)/lib
```

### Option 3: Use Docker

First build the container:
```
docker build -t surge:latest .
```

Then run it by passing arguments directly (the image uses an ENTRYPOINT):
```
docker run --rm surge:latest -u C10H16O
docker run --rm surge:latest -S C6H6
```

To save output to a file on the host:
```
docker run --rm surge:latest -S C10H16O > results.smi
```

## Misc
Surge was developed by [Brendan McKay](http://users.cecs.anu.edu.au/~bdm) with the help of [Christoph Steinbeck](https://github.com/steinbeck) and [Mehmet Aziz Yirik](https://github.com/mehmetazizyirik).

