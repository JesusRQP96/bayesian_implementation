# Analisis bayesiano
* Analisis espacial de la BD  Pennsilvania : un enfoque bayesiano 
* Analisis multinivel de brecha de genero : Un enfoque espacial

# Analisis multinivel de brecha de genero
El modelo implementado en este estudio fue uno de corte bayesiano puesto que dicho enfoque tendria principalmente 2 ventajas :
1. Modelar el ingreso como variable latente (beneficioso por como los datos de ingresos fueron registrados)
2. Facilitar la modelacion jerarquica, más aun cuando el objetivo del estudio era evaluar efectos a nivel departamental

El modelo implementado fue finalmente el siguiente,el presente repositorio contiene los codigos base para implementar un modelo similar en *stan* aunque de menor complejidad 
![equation](http://www.sciweavers.org/upload/Tex2Img_1603147732/render.png)

## Principales resultados 
Se obtuvieron resultados para el ingreso promedio y brecha salarial, como los siguientes  :

![alt text](https://github.com/JesusRQP96/bayesian_implementation/blob/master/multinivel-med-avg-wage-figure.png)
![alt text](https://github.com/JesusRQP96/bayesian_implementation/blob/master/multinivel-med-wage-gap-figure.png)




