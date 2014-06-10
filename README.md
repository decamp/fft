### Fast Fourier Transform Library:
There are many like it, but this one is mine. 

This is a very simple library that's fairly well optimized. 
It contains pure java implementations of:
- 1D Fast Fourier Transform
- 2D Fast Fourier Transform
- 1D Fast Cosine Transform  
- 2D Fast Cosine Transform

1D transforms only operate on vectors where the length is a power-of-two.
2D transforms only operate on square matrices where the size is a power-of-two.


### Build:
$ ant


### Runtime:
After building, add all jars in **target** directory to your project.


### Dependencies:
None.

---
Author: Philip DeCamp
