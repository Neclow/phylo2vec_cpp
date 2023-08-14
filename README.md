# Phylo2Vec: a vector representation for binary trees

## Installation
Prerequisites:
 * C++14
 * GoogleTest 1.11.0: ```sudo apt-get install libgtest-dev```
 * clang-format: ```sudo apt install clang-format```
 * cmake 3.22.1: ```sudo apt-get install cmake```
 * cxxopts 3.1.1:
    ```
    git clone https://github.com/jarro2783/cxxopts.git
    cd cxxopts
    mkdir build
    cd build
    cmake ..
    make
    sudo make install
    ```

Compile and build from scratch:
```
mkdir build
cd build
cmake ..
make
```

Or use the executable. To get information about all arguments:
```
./phylo2vec -h
```

Arguments:
```
  -h, --help            Show help
      --toNewick arg    Convert to Newick format. Example input: 0 1 4
      --toVector arg    Convert to integer vector. Example input:
                        "(((2,1)4,0)5,3)6;"
      --num_leaves arg  Number of leaves (optional)
```

Example usage of toNewick:
```
phylo2vec --toNewick 0 1 4
```

Example usage of toVector:
```
phylo2vec --toVector "(((2,1)4,0)5,3)6;" --num_leaves 4
```

## TODO list:
* toVector: Accept a newick with string-formatted nodes (not just ints) and output the vector with a int->string mapping

## Python version:
* https://github.com/Neclow/phylo2vec
* A minimalistic demo is available on Colab: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/10ZENm-wgWiRFa4ABY8piGDY_QoJyZ30X?usp=sharing)

## Citation and other work
* Link to the preprint: https://arxiv.org/abs/2304.12693
* Related work: https://github.com/Neclow/GradME = phylo2vec + minimum evolution + gradient descent
