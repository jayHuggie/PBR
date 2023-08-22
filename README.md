# A Physically Based Renderer: Course Project of UCSD CSE 168

For the course project for ‘UCSD CSE 168 Computer Graphics II: Rendering’, we constructed a physically-based renderer.

Detailed feature overviews and conceptual breakdowns can be found in our assignment specification:


# Features:

* Parses a scene in Mitsuba XML file format.
* Path Tracing
* Area Light (Monte Carlo Integration)
* Multiple Importance Sampling (we adopted the one-sample variant of MIS: instead of deterministically shooting rays for both lights and BRDFs and weighing them, we randomly choose one and combine the distribution)
* Implements a BVH acceleration structure.
* Support for spheres and triangle meshes.
* UV Texture Mapping
* Support for diffuse, mirror, plastic materials and Phong, Blinn-Phong, Blinn Phong microfacet BRDF
* Anti-aliasing, Russian Roulette termination.
* Exports the scene into a EXR file.

# Build

1. Clone the repository.
2. Navigate to the 'build' directory.
3. Use CMake to build (this instruction assumes you're on a Unix-like system).
4. Use the 'make' command to compile the code. Note: This project was tested on an M1 MacBook Air.


```
git clone https://github.com/jayHuggie/PBR
cd PBR/build
cmake ..
make
```

# Running the Code

After compiling, use the command:
```
./torrey -hw 4_3 ../scenes/veach_mi/mi.xml
```

* Regarding the '-hw' argument: It accepts values from 1_1 to 4_3, representing different functions detailed in the assignment documents. For the most comprehensive experience, use '-hw 4_3'.  <br />
* '../scenes/veach_mi/mi.xml' after -hw 4_3 is just an example. Additional scenes can also be rendered. The xml file can be located in various folders within the ‘scenes’ folder like:
```
./torrey ..scenes/cbox/cbox.xml
```
* The resolution and sample size can be adjusted by changing the values of the xml files.
