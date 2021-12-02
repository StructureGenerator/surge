Surge.exe is a command line application that you need to run from the Windows command line by typing "surge.exe" and hit return.
Please make sure that the zlib1.dll is in the same folder as surge.exe.
This dll enables surge to store the different output formats in a compressed manner by using the -z command line option.
surge.exe -help will list all available options.
For example "surge.exe -u C10H16O" will generate all 452458 isomers of this formula but only report the final count and not output the structures.
"Surge.exe -S C10H16O" will do the same but also output the structures in SMILES format.
"Surge.exe -Sz C10H16O > C10H16O.smi.gz" will compress the SMILES and write the output into a file. 
