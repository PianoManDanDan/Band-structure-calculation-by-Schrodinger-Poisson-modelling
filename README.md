# Band Structure Calculation by Schrodinger-Poisson modelling

This research project was done as part of my masters year at Cardiff University.

### Dependecies:

- Python 3 (v3.7.0 was used in development) - While efforts have been made to ensure this is compatible with Python 2, it has not been tested so may not work as intended.
- [Numpy Library](http://www.numpy.org/) (v1.15.2 was used in development)
- [Matplotlib Library](https://matplotlib.org/) (v1.1.0 was used in development)
- [Scipy Library](https://www.scipy.org/) (v3.0.0 was used in development)

### How to use

To run the Schrodinger-Poisson solver:
1. Run `main.py`.
2. When the window appears, type in a temperature (in Kelvin) in the box at the top of the left hand side.
3. Select a material from the dropdown menu. All of the values for the material should appear in the boxes along the rest of the row. 
* Note: If you have selected Al<sub>x</sub>Ga<sub>(1-x)</sub>As, you will also need to enter an x value (between 0 and 1) before these values will appear.
* If you have selected 'Custom' you will need to manually enter each of the values.
4. Enter a layer thickness (in nm).
5. Click on the '+' button to add another layer. A minimum of 3 layers need to be added before you can calculate the material potential. You cannot add further layers until all previous layers have been filled in. You can add up to a maximum of 25 layers.
6. Click on the 'Calculate Potential' button to calculate the potential and plot it on the figure on the right hand side.
7. Click on the 'GO!' button to calculate the Schrodinger-Poisson solutions.

You can adjust the number of wavefunctions to find, as well as the number of points used by changing the Settings in the 'Settings' tab. The default settings are to find 5 wavefunctions for each potential and to use 1000 points.

