# Band Structure Calculation by Schrodinger-Poisson modelling

This research project was done as part of my masters year at Cardiff University.

### Dependecies:

- Python 3.6 or later
- [Numpy Library](http://www.numpy.org/) (v1.15.2 was used in development)
- [Matplotlib Library](https://matplotlib.org/) (v1.1.0 was used in development)
- [Scipy Library](https://www.scipy.org/) (v3.0.0 was used in development)

### Setup

To run the Schrodinger-Poisson solver, run `main.py`, select a material and potential profile, and click on the 'GO!' button.

### Choosing your own material

There are several built in preset materials, which you can select from the drop down menu. Alternatively, you can select your own material properties by  selecting the appropriate option from the dropdown menu. When creating your own material profile, the variables should be stored in a `.json` file in the style as given below:

```json
{
	"me": "float: Effective mass of electron (kg)",
	"mhh": "float: Effective mass of heavy hole (kg) ",
	"mlh": "float: Effective mass of light hole (kg)",
	"Eg": "float: Band Gap of material (eV)",
	"dielectric_constant": "float: dielectric constant of material",
}
```

The variable names `"me"`, `"mhh"`, `"mlh"` and `"Eg"` must be surrounded by quotation marks (**"**), followed by a colon (**:**) and then a floating point number (without quotation marks around it). Each line must terminate with a comma (**,**).

The example above can be used as a template, and can also be found in the `./materials/template.json` file.

### Creating a potential

To choose a potential, the spatial values and corresponding potential values should be stored in a `.csv` file **with spatial values in the first column, and potential values in the second**. 

For example, a simple, infinite potential well between -10nm and 10nm with 1000 points can be created using this simple python script:

```python
import numpy as np 

x = np.linspace(-10e-9, 10e-9, 1000) # metres
V = np.zeros_like(x) # Joules

potential_profile = np.array([x, V])
potential_profile = np.transpose(potential_profile) # Needed to make column values

np.savetxt('potential.csv', potential_profile, delimiter=',')

```

**The number of spatial points must be the same as the number of points in the potential.**

