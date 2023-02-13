# multicore_project

program for image filtering. lets you blur, sharpen, find edges and leak memory, both in a fast cuda accelerated version and in a slow cpu based version. the fast version resides in multicore_project/PROJECT/filters and the slow version is located in multicore_project/PROJECT/filtersSlo. obviously the meat of the source code are the files with the .cu extension.

I have compiled them using 'make SMS="75"' with the directory set to /filters or /filtersSlo on linux, on one of the computers located in room 204 of the faculty of physics of AGH. prior to this a 'source setcuda' command was input while in the /multicore_project directory.

good luck to anyone attempting to do this anywhere else and especially on windows. from what I found compiling either of the .cu files with nvcc works just as well for whatever reason even without any magic arguments and options. filtersSlo.cu, in spite of its extension, is a pure C++ program and may be compiled using g++ or other C++ compilers and ran without cuda hardware. I have stopped doing this in the latest testing phase so that the programs could be directly compared with a clear conscience.

the folders of programs contain the source codes and compiled programs as were used for composing the report, which is located at multicore_project/multicore_project.pdf and probably is a reasonably good read if you find these programs of any interest.

in the program folders are other readme files which dive into more intricate (less understandable) details of the project.

there is a sample input folder containing some samples and a sample output folder as well. also included is the cuda environment which does many important things.

thanks for clicking here or maybe even reading
