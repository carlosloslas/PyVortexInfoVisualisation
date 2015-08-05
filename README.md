## Vortex Information Visualisation with Python
#### Repository containing the code from my part III Ship Science disertation.

### Introduction
The interaction between a fluid and a structure is a phenomenon very relevant to the field of aero-hydrodynamics. When the intection is modeled in detail there is presence of vortical structures. There is a need for informative visualisation of these vortical structures in order to fully understand the fluid-structure interaction.

*Jeong and Hussain* [1] proposed a strong mathematical definition for the core of a vortex. This project took this definition and developed a range of algorithms to beter visualise the vortices within the velocity field.

In order to simplify the problem the fields where self-generated. But for the algorithm to have value, it was developed form a field-agnostic point, that is that it would work for any field with a number of vortices in it. 

The wake of a fish was chosen as a flied to generate. Which is defined as a number of vortex rings in alternating directions. The reason begind choosing this wake was due to the clarity and simplicity of its vortical structures.

The results obtained show that there is an iprovement in the visualisation of the vortex rings. The starting point does not provide any information other that the magnitude and direction of the velocity field. When the definition for the vortex core [1] is aplied the rings are revealed, providing the shape of the vortical structures. During this project an algorithm capable of obtaining the circulation was developed, which added meaning to the initial shape provided by the mathematical definition [1].

As a breif conclusion it can be said that the project contributed to improve the visualisation techniques of vortical flow fields. Allthough the project didn't achieve visualisation of the vorticity of the vortex rings, this is suggested in the future work section of the full write-up of the project, together with code and algorithm improvements. The work initially developed could lead to further research into an efficient way of visualising the vorticity of any given flow field where vortices are present.

### Algorithm elements
1. Velocity field generation. Functions which calculate and display the velocity field induced by one or more vortex rings.

2. VTK image writing. Functions that take the generated velocity field and translate it to a binary VTK image using the package 'tvtk' package.

3. Image reading. Using the 'vtk' package to read the binary image where the field is stored.

4. Vortical field information extraction and visualisation. Combining the lambda 2 vortex core definition, linear algebra, and vortex properties, to firstly locate the vortex core and subsequently extract the circulation from the field. 

### Results obtained

#### Vortex rings
Self-generated velocity field due to a single vortex ring.
![3dvort](https://cloud.githubusercontent.com/assets/10100481/9084759/9c821fe0-3b6e-11e5-98d6-9aec0611dea5.png)
Aplication of the Lambda 2 definition of a vortex core.
![3dl](https://cloud.githubusercontent.com/assets/10100481/9084764/a2ab814a-3b6e-11e5-87b7-494747a63923.png)
Use of lienear algebra and vortex properties to obtain and visualise the circulation of the vortex ring.
![ringc](https://cloud.githubusercontent.com/assets/10100481/9084767/a8001de0-3b6e-11e5-9ac7-c9a1c68610ac.png)

#### Vortex ring chains
Self-generated velocity field.
![3dvchain](https://cloud.githubusercontent.com/assets/10100481/9027607/cf689992-3953-11e5-8f24-4a60161b2c20.png)
Lambda 2 definition of the vortex core.
![chainl](https://cloud.githubusercontent.com/assets/10100481/9027612/01a8de30-3954-11e5-92a0-efbd0d02b373.png)
Circulation of the vortex rings.
![chainc](https://cloud.githubusercontent.com/assets/10100481/9027614/0b578936-3954-11e5-9b86-479a22d94fb7.png)

### How to use it?

1. I developed this project using the Anaconda distribution of Python. I can't guaranty that it will work if you dont use it. http://continuum.io/downloads.
2. Firstly clone the repository into your system. 
`git clone https://github.com/carlosloslas/PyVortexInfoVisualisation.git`

3. Open Spyder, the editor I used for the developement of the project. In your terminal type: 
`spyder`
4. Setup your spyder working directory to where you cloned the repository.
5. Open the ``` visual.py ``` file.

###Bibliography
[1] Jinhee Jeong and Fazle Hussain, On the identification of a vortex, Journal of Fluid Mechanics, Vol. 285 (1995), pp. 69-94.

### Acnowledgements

Thanks to my supervisor Dr. Gabriel Weymouth for all his support. Without all his good suggestions and meetings this proyect would have never reached this level. 
