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

Example usage of toVector:
```
./phylo2vec --with_mapping --toVector "(((((((tip_0:1.44,tip_1:1.44)8042:0.46,(tip_2:1.5,tip_3:1.5)8043:0.4)8044:0.3,(tip_4:1.51,tip_5:1.51)8045:0.69)8046:0.4,tip_6:2.6)8047:1.05,tip_7:3.65)8048:0.5,(((tip_8:0.72,tip_9:0.72)8049:0.28,tip_10:1)8050:1.56,tip_11:2.56)8051:1.59)8052:1.96,tip_12:6.11)8053:0;"
```

## Python version:
* https://github.com/Neclow/phylo2vec written with [Matthew Penn](https://www.stats.ox.ac.uk/people/matthew-penn) and [Samir Bhatt](https://publichealth.ku.dk/about-the-department/section-epidemiology/?pure=en/persons/707469)
* A minimalistic demo is available on Colab: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/10ZENm-wgWiRFa4ABY8piGDY_QoJyZ30X?usp=sharing)

## Citation and other work
* Link to the preprint: https://arxiv.org/abs/2304.12693
* Related work: https://github.com/Neclow/GradME = phylo2vec + minimum evolution + gradient descent
* The file ```100trees.txt``` is an anonmyised subset of [TimeTree](http://www.timetree.org/) trees adapted from: https://datadryad.org/stash/dataset/doi:10.5061/dryad.2fqz612vg. Kudos to [Mark Khurana](https://forskning.ku.dk/soeg/result/?pure=da%2Fpersons%2Fmark-poulsen-khurana(171ece7e-9567-4d48-8cf9-959b57de57c8).html) for the dataset

## Note
* At least half of this repo was coded on a Samsung A40 using [JDoodle](jdoodle.com), that should probably deserve a star :)