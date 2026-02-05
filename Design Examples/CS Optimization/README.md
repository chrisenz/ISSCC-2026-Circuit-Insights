# Optimization of a Common-source Stage in Open-loop

![CS OL Amplifier](/Design%20Examples/CS%20Optimization/Figures/CS_OL_schematic.png)

This [Quarto notebook](/Design%20Examples/CS%20Optimization/CS_optimization.pdf) shows an example of how to design a CS stage for minimum power consumption in open-loop configuration using the inversion coefficient. The circuit is designed for the IHP 130nm BiCMOS process using the sEKV parameters extracted in the companion notebook. A first example is showing the optimization of the CS stage for a given gain-bandwidth product and a given transistor length. The second example is taking advantage of the additional degree of freedom on the transistor length to simultaneously optimize for a given DC gain and gain-bandwidth product.
