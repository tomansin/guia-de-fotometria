import numpy as np
import os
import subprocess
import sys
from IPython.display import clear_output

from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from ccdproc import ImageFileCollection

def runCommand(cmd, timeout=None):
    ''' 
    La siguiente función sirve para correr comandos de bash desde el notebook 
    '''
    
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = ''
    for line in p.stdout:
        line = line.decode(encoding=sys.stdout.encoding,
                    errors='replace' if (sys.version_info) < (3, 5)
                    else 'backslashreplace').rstrip()
        print(line)
        output += line

    retval = p.wait(timeout)
    #log.debug('retval=%d' % retval)

    return (retval)

def Astrometria(imagenes, dir_salida, escala, n=20,forzar_escritura=False):
    '''
    Esta función agrega información asrtometríca a una serie de imágenes,
    y las guarda en una carpeta. Ademas devuelve una tabla con información
    sobre las imágenes.
    
    Parámetros 
    ----------
    imagenes: `~ccdproc.ImageFileCollection`
        Colección de las imágenes que se quiere procesar
    dir_salida: str
        Nombre del directorio de salida
    escalas: list
        Lista de las escalas [scala_low, scala_up] para dicho telescopio
    n: int    
        Número de fuentes mas brillantes utilizado para estimar el fwhm
    '''
    
    # Armo tabla de estadística    
    stats_table=Table(names=('fname','fwhm','fwhm_std','img_mean','img_std'),
                  dtype=('a15','f','f','f','f'))
    
    
    # Armo el encabezado y el pie de la nueva sección del header
    encabezado = (('COMMENT' ,'Original key: "END"'), 
              ('COMMENT','--Start of Astrometry.net WCS solution--'), 
              ('COMMENT'))
    pie = (('COMMENT'), ('COMMENT','--End of Astrometry.net WCS--'), 
           ('COMMENT')) 
    
    if forzar_escritura:
        print('Las nuevas imágenes seran guardadas en {}'.format(dir_salida))
    else:
        if os.path.exists(dir_salida):
            raise FileExistsError("Para escribir en una cerpeta existente "
                                  "usar el parametro 'forzar_escritura=True'")
        else:
            os.mkdir(dir_salida) 
            print('Se creo el directorio {}, que contendra las nuevas imágenes'.format(dir_salida))

    # Astrometry necesita valores limites de escala
    if not 0 < escala <= 1:
        raise ValueError("'escala' debe estar en el rango (0,1)")

    escala_up = str(np.ceil(escala*10)/10)
    escala_low = str(np.floor(escala*10)/10)
    
    # Itero por cada una de las imágenes en el ImageFileCollection
    for obj, fname in imagenes.hdus(imagetyp='object',return_fname=True):
        clear_output(wait=True)
        print('Buscando fuentes en imagen:', fname)

        # guardo header y data
        hder = obj.header
        data = obj.data
        
        # modifico el header para que sea fits standar
        del hder['EPOCH']
        del hder['EQUINOX']
        
        # por alguna razon el DATE-OBS no coincide con MJD-OBS para HSH
        hder['DATE-OBS']= Time(hder['MJD-OBS'], format='mjd').isot
        
        # extraigo valores de RA y DEC del header para astrometry
        img_coo = SkyCoord(hder['RA2000'], hder['DEC2000'],
                       unit=(u.hourangle, u.deg))
        ra, dec = str(img_coo.ra.degree), str(img_coo.dec.degree)
    
        # extraigo tamaño del ccd para astrometry
        width = str(hder['naxis1'])
        height = str(hder['naxis2'])
        
        # corremos sextractor -> genera test.cat (archivo binario)
        runCommand('sextractor '+imagenes.location+fname)
        coords = 'test.cat'
        t = Table.read(coords)
        rad_mean = np.mean(np.sort(t['FLUX_RADIUS'])[::-1][:n])
        rad_std = np.std(np.sort(t['FLUX_RADIUS'])[::-1][:n])
        img_mean,_,img_std = sigma_clipped_stats(data)
        stats_table.add_row([fname,rad_mean*2,rad_std*2,img_mean,img_std])
        
         
        clear_output(wait=True)
        print('Resolviendo astrometría a imagen:', fname)
    
        # corremos astrometry -> genera .wcs con info astrométrica
        runCommand('solve-field --overwrite --temp-axy \
                    --index-xyls none --solved none --match none --rdls none \
                    --corr none --no-plots --ra '+ra+' --dec '+dec+' \
                    --radius 1 --scale-units arcsecperpix --scale-low '+escala_low+' \
                    --scale-high '+escala_up+' --x-column X_IMAGE --y-column Y_IMAGE \
                    --sort-ascending --sort-column MAG_AUTO --width '+width+' --height '+height+' \
                    --wcs obj.wcs '+coords)
    
        clear_output(wait=True)
        print('Nombre de imagen:', fname)

        # armo las nuevas imágenes con la información astrometrica
        # primero abro obj.wcs y obtengo el header como una lista
        with fits.open('obj.wcs') as whdu:
            whdr = whdu[0].header
            wlist = list(dict.fromkeys(whdr))
        # agrego informacón astrométrica al header de la imagen
        # junto con encabezado y pie
        wlist.remove('DATE') # borro keyword repetido
        for enc in encabezado:
            hder.append(enc, end=True)
    
        for key in wlist[4:]:
            if key == 'HISTORY':
                for his in whdr[key]:
                    hder.append(('HISTORY',his), end=True)
            elif key == 'COMMENT':
                for com in whdr[key]:
                    hder.append(('COMMENT',com), end=True)
            else:
                hder.append((key,whdr[key],whdr.comments[key]), end=True)

        for p in pie:
            hder.append(p, end=True)

        # genero nueva imagen
        fits.writeto(dir_salida+'/'+fname[:-4]+'w.fits', data, hder, overwrite=True)

        clear_output(wait=True)

    clear_output(wait=True)
    print('Las nuevas imágenes con astrometría fueron guardadas en la carpeta {}'.format(dir_salida))
    
    # Guardo la tabla de la estadística
    stats_table.write('images_prop.dat',overwrite=True,format='ascii')
    
    # Borro archivos que no me sirven
    os.remove('test.cat')
    os.remove('obj.wcs')