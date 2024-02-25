Moody installation

Moody itself is released in pre-compiled form, however the API coupling code
is released open source. Therefore (as of v.2.0.1), only the binaries are packed
in .tar.gz files for each operating system (OSX, Linux and Windows). The com
pressed files includes the etc, bin, lib and include folders. The tutorials, the
manual and the API are all accessed directly as sources when cloning the reposi
tory.

  1. Clone or download the repository from github.com/johannep/moodyAPI.
  2. Select your installation directory and unpack the appropriate tar-file bina
     ries IN THE SAME LOCATION. Technically, this is only relevant for all
     tutorials and matlab-processing paths to work as intended.
  3. In Linux or OSX environments, source the environment variables. For this,
     use etc/bashrcMoody. NB: There is no Windows PATH-setting script
     provided at this time.

Ex: type these lines in the command window ù8Ubuntu or centOS)

 cd moodyAPI;
 
 tar-xzvf moody-Linux.tar.gz;
 
 source etc/bashrcMoody

Once the previous steps are completed, the moody postProcessing utility can be used. moodyPost.x can be used to create a smaller
 set of output data for post-processing. It can also be used to generate VTK-files of
 the cable lines and their tension force magnitude for visualisation in e.g. Kitware’s
 Paraview. 

 The postprocessing utility can be used to save -vtk files as well as files containing for each integration time-step the values of tension, strain, position, and velocity.
