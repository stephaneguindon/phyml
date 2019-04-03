
### 1. You downloaded a stable release  of PhyML from [here](https://github.com/stephaneguindon/phyml/releases) 

To install any program that is part of the PhyML package, type the following commands:

```bash
./configure --enable-XXXX;
make;
```
where XXXX is phyml or phyrex or phytime.

To compile a Windows executable, install MinGW and run:

```bash
./configure --enable-win --enable-XXXX;
make;
```

To install the MPI version of PhyML, type the following commands:

```bash
autoreconf -i;
./configure --enable-mpi --enable-phyml;
make;
```

If you are using a Mac computer, you will need to install the package pkg-config. The
following command should set you up (provided Homebrew is installed on your Mac...):
brew install pkg-config; Next, typing `./configure --enable-XXXX; make;` should generate
the phyml binary in the `src/` directory.



### 2. You have cloned PhyML from GitHub

To install any program that is part of the PhyML package, type the following command:

```bash
sh ./autogen.sh;
```

If you are using a Mac computer or running a Unix-like operating system, you will need 
to install the packages autoconf automake and pkg-config. On a Mac, the following command 
should set you up (provided Homebrew is installed on your Mac...): brew install pkg-config
autoconf automake;

Next, typing `./configure --enable-XXXX; make;` should generate the phyml binary in the src/
directory, where XXXX is phyml or phyrex or phytime.


