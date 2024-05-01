# chem279_final_project

Step 1:
Add unit test checking that python was correctly getting points from the c++ code
    1a: Determine point generation method
    1b: Write code for point evaluation

Step 2: 
Make a few example python plots with MatPlotLib to practice & demonstrate plotting techniques with matplotlib

Step 3:
Refactor previous class code to my purposes and preferences
    - Side note: personally I'm pretty big on no abbreviated variable names when I'm authoring something I lose confidence in this frequently haha

Step 4: 
Iterate through many different plotting techniques 
Top 2 choices currently: 
    marching_cubes algorithm
    scatter (but panning has a bug)

Step 5: 
Plot Molecular orbitals
    5a: Decide visual representation of scalar MO coefficients

## Code conventions
Molecule/Atom symbols will be capitalized: Class names will be capitalized, pretty much everything else will be lower case and underscore separated.

## References

Marching Cubes Algorithm for Isosurface plotting:
https://scikit-image.org/docs/stable/auto_examples/edges/plot_marching_cubes.html

3d plotting in matplotlib -
https://www.youtube.com/watch?v=fAztJg9oi7s
https://www.youtube.com/watch?v=PnwpoCDA5IM

## Problems encountered

Open bug exists in matplotlib 3d scatter plot opacity:
https://stackoverflow.com/questions/71904575/matplotlib-3d-scatter-plot-alpha-varies-when-viewing-different-angles
references https://github.com/matplotlib/matplotlib/issues/22861

Challenging to keep grid transformations straight between c++ and python grids/arma::cubes/np arrays, linspaces for subsequent use in matplotlib, etc

# Usage
To compile the code, enter the command `make all` under this directory.

```
./HW5 <filename>
```
`<filename>`: The path to the input file containing information of the molecule.

