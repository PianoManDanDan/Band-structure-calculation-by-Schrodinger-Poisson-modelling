# Band Structure Calculation by Schrodinger-Poisson modelling

This research project was done as part of my masters year at Cardiff University.

### Dependecies:

- Python 3.6 or later
- [Numpy Library](http://www.numpy.org/) (v1.15.2 was used in development)
- [Matplotlib Library](https://matplotlib.org/) (v1.1.0 was used in development)
- [Scipy Library](https://www.scipy.org/) (v3.0.0 was used in development)

### Setup

To run the Schrodinger-Poisson solver, run `main.py`. When the window appears, select a material from the drop down box and type the thickness of that growth layer and the temperature of the material. For AlGaAs, also select the Aluminium percentage (given as a float between 0 and 1). Click on 'Calculate Potential' to calculate the potential and give a visualisation of the chosen potential. Press 'GO!' to run the Schrodinger-Poisson solver.

To add more growth layers, click on the '+' button.

### Choosing your own material

There are several built in preset materials, which you can select from the drop down menu. Alternatively, you can select your own material properties by  selecting the 'Custom material' option from the dropdown menu. When creating your own material profile, you must enter all values into the appropriate boxes.

### Creating a potential

When you select a material, a thickness and a temperature, you can click on 'Calculate Potential' to calculate the potential profile of your growth material. There should be a minimum of 3 materials chosen - two substrate layers and a growth layer - however, several growth layers can be added on to create the potential.

