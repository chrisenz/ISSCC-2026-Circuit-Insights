# ISSCC-2026-Circuit-Insights
This repository gives additional material to the [ISSCC 206 Circuit Insights](https://www.isscc.org/insights) course which detailed program is given in the [ISSCC 2066 Circuit Insights Flyer](https://www.isscc.org/s/2025-12-02-ISSCC2026-Circuit-Insights-Flyer.pdf)

In this repository you will find the detailed [sEKV parameters extraction procedure](/sEKV%20Parameter%20Extraction/README.md) for the [open source IHP SG13G2 130nm BiCMOS technology](https://github.com/IHP-GmbH/IHP-Open-PDK). You also find an example showing how to design a [common-source stage](/Design%20Examples/CS%20Optimization/README.md) for minimum power consumption using the inversion coefficient for the same IHP 130nm technology.

If you are only interested to learn how the sEKV parameters are actually extracted you can just have a look at the pdf files. If you want to perform the extraction yourslef you can run any of the Quarto notebboks. For doing this you need to have the same [installation](/Installation.md) as me.

If you want to run the extraction on different transistor sizes, you need to generate the required data using for example ngspice. For doing this you need to install ngspice following the [ngspice installation instructions](ngspice_installation.md).

Finally, if you want to extract the sEKV parameters on a different technology you can reuse the Quarto notebooks, but you then need to do have the PDK of the target technology.

I hope you will find this repository useful! Any feedback is welcome. Send me a mail at christian.enz@epfl.ch.

Enjoy!

Christian Enz
