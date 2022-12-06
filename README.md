# SAM reader

Used to extract simple data from a .sam file.

## Description

A Python script that extracts simple data from a .sam file. It needs a .sam file (see [samformat.info](https://www.samformat.info/) for more info) in parameter. When running the script, it will ask the user to input a certain threshold for filtering reads with a certain mapping quality. It will generate a .txt file that sums up simple data.

Input:
* a single .sam file

Output:
* a .txt file containing:
    * the total number of reads
    * the number of different reads
    * the number of paired reads
    * the number of reads mapped in proper pair
    * the number of reads matched at 100% (CIGAR)
    * The number of reads with a mapping quality higher than a certain threshold (threshold can be set manually)

## Getting Started

### Dependencies

* Linux-derived OS
* Python3

### Executing program

```
./sam_reader.py sam_file.sam
```

## Authors

Elian Strozyk <elian.strozyk@etu.umontpellier.fr> at Université de Montpellier

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>."

## Version History

* 0.0.1 (2022-12-06)
    * Initial Release

## Acknowledgments

* Anna-Sophie Fiston-Lavier
* Arnaud Soulier and Allyson Moureaux
