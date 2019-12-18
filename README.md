# Guia de fotometría diferencial
Esta guia consta de una serie de notebooks para realizar fotometría diferencial, utilizando **SEXTRACTOR**, **Astrometry.net** y las librerías de Python: **Numpy**, **Matplorlib**, **Pandas**, **Astropy** y **Photutils**.

## Como usar

Partiendo de una serie temporal de imágenes astronómicas con formato FITS

* **ASTROMETRÍA:** Las imágenes deben tener solución astrométrica. Para eso se procede con el notebook **Astrometria.ipynb**. Es importante que dichas imágenes esten en formato FITS estandar, que tenga los keywords 'RA2000' y 'DEC2000', y que el directorio 'sextractor' se encuentre en el directorio de trabajo.

* **FOTOMETRÍA:** Para realizar la fotometría se usa el notebook **Fotometria.ipynb**. Creará un directorio con la fotometría para cada imágen  con los detalles del procedimiento. Por otro lado devolverá una tabla con las coordenadas celestes de las fuentes, y otra tabla con la fotometría como serie temporal.

* **CURVAS DE LUZ:** La fotometría diferencial se obtiene luego de una serie de pasos detallados en el notebook **Analsis.ipynb**.

## Instalación de dependencias

La instalación de **Astrometry.net**, instala **SEXTRACTOR**.
Para instalar **Astrometry.net** en distribuciones basadas en Debian:

```bash
$ apt install astrometry.net
```
De otra manera puede instalarlo desde la pagina:
http://astrometry.net/use.html

Para instalar las librerías de Python en conda:
 
```bash  
$ conda install numpy    
$ conda install matplotlib  
$ conda install pandas  
$ conda install astropy
$ conda install photutils -c astropy
$ conda install ccdproc -c astropy
$ conda install reproject -c astropy
```

