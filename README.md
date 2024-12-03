# Genomics Parallel Edit Distance Project

## Repo Setup

On your machine you should have at least
- g++/gcc 11
- CMake 3.22.1
- Ninja 1.10.1
- Git 2.34.1

After cloning the repo, run 
```
git submodule update --init --recursive
```

To setup VCPKG
```
cd vcpkg && ./bootstrap-vcpkg.sh
export VCPKG_ROOT=vcpkg
export PATH=$VCPKG_ROOT:$PATH
```

To build with cmake
```
cmake --preset=default
cmake --build build
```

## Scripts

input file is two strings on two lines with newlines seperating the two lines
```
<STRING_1>
<STRING_2>
```

Penalty Function CSV is CSV with three columns: original character, changed character, penalty. Gaps are represented as -
```
A,A,0
A,C,4
A,G,2
A,T,4
A,-,8
...
```

### Hirschberg

```
./build/hirschberg <input_file> <penalty_function_csv>
```

### Global/Local Alignment
For Affine Global Alignment
```
./build/alignment <input_file> affine <iterating_mode> <penalty_function_csv> <gap_open> <gap_extension>
```
For Local/Global Alignment
```
./build/alignment <input_file> <global|local> <iterating_mode> <penalty_function_csv>
```

### Block Based Global/Local Alignment

For Affine
```
./build/alignment <input_file> affine 3 <internal_iterating_mode> <block_iterating_mode> <block_size> <penalty_function_csv> <gap_open> <gap_extension>
```
For Local/Global Alignment
```
./build/alignment <input_file> <global|local> 3 <internal_iterating_mode> <block_iterating_mode> <block_size> <penalty_function_csv>
```
### Four Russians Based Edit Distance
```
./build/hirschberg <input_file> <penalty_function_csv>
```
