# sallow
A heuristic algorithm for finding [treedepth decompositions](https://en.wikipedia.org/wiki/Tree-depth).

This is a submition for the [PACE 2020](https://pacechallenge.org/2020/td/) challenge.
Sallows (such as ([Salix caprea](https://en.wikipedia.org/wiki/Salix_caprea)) are a kind of willow shrubs/trees. 

![A sallow tree](https://upload.wikimedia.org/wikipedia/commons/thumb/c/c2/20170421Salix_caprea1.jpg/640px-20170421Salix_caprea1.jpg)

(Photo by [AnRo0002](https://commons.wikimedia.org/wiki/File:20170421Salix_caprea1.jpg))

## Git
This repo uses submodules:
* clone with `git clone --recurse-submodules git@github.com:marcinwrochna/sallow.git`
* update all submodules with `git submodule update --remote --rebase` (this fetches to `master` and tries to rebase your existing changes within each submodule)

## Building
Standard CMake:
* `mkdir build && cd build` (or wherever you want the build files to be, instead of `build/`)
* `cmake --config Release ..` (or whatever the path to the root CMakeLists.txt is, instead of `..`)
* `make`
* Or just open the directory in [MSVS Code](https://code.visualstudio.com/) with the [CMake extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cmake-tools) installed, and allow it to configure itself.

## Running
* `sallow input.gr`
* Interrupt (ctrl+C or SIGINT) to stop and print the best decomposition we have.
* If it hangs (e.g. you ran it in debug mode on extremely large graphs), you need to `pkill -9 sallow`.
* Input and output format as specified [here](https://pacechallenge.org/2020/td/) (DIMACS-like).
