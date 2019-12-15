import numpy as np
import pandas as pd
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
from astropy.wcs import WCS

from ccdproc import ImageFileCollection
from photutils import aperture_photometry,SkyCircularAperture, SkyCircularAnnulus, DAOStarFinder


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
     
    
    
def Fotometria(imagenes, aperturas, anillos, dir_tablas='Tablas', forzar_escritura=False):
    '''
    Esta función realiza la fotometría a una serie de imágenes, dadas las aperturas,
    en coordenadas celestes.
    
    Parámetros
    ----------
    imagenes: `~ccdproc.ImageFileCollection`
        Colección de las imágenes que se quiere procesar
    aperturas: `~photutils.aperture.circle.SkyCircularAperture`
        Aperturas en coordenadas celestes.
    anillos: `~photutils.aperture.circle.SkyCircularAnnulus`
        Anillos para calcular el cielo. Tambien en coordenadas celestes
    dir_tablas: str, opcional
        Directorio donde se guardaran las tablas de fotometría de cada imágen
    forzar_escritura: bool, opcional
        Forzar toda la escritura si los archivos existen. Tener cuidado, se puede
        perder datos
        
    Devuelve
    --------
    fot_tabla: `~astropy.table.Table`
        Tabla de fotometría
    '''
    
       
    if forzar_escritura:
        print('Las nuevas tablas quedaron guardadas en el directorio {}'.format(dir_tablas))
    else:
        if (os.path.exists(dir_tablas)) & (os.path.exists('fotometria.csv')):
            raise FileExistsError("Para sobreescribir usar el parametro "
                                  "'forzar_escritura=True'")
        else:
            os.mkdir(dir_tablas) 
            print("Se creo el directorio {}, que contendra las "
                  "nuevas tablas".format(dir_tablas))
            
    ###########################################################################
    
    zmag=25 #defino magnitud arbitraria tal como IRAF
    
    # Función para calcular el error en la fotometría. Igual que IRAF
    def error(flux, gain, area, stdev, nsky):
        return np.sqrt((flux/gain) + area*stdev**2 + area**2*stdev**2/nsky)
    
    # Armo la tabla de fotometría
    colnames=['FNAME', 'MJD','DATE']
    for i in range(0,len(aperturas)):
        mag = 's'+str(i+1)
        merr = 'merr'+str(i+1)
        colnames.append(mag)
        colnames.append(merr)
    fot_tabla = Table(names=colnames)
    fot_tabla['FNAME'].dtype = 'a15'
    fot_tabla['DATE'].dtype = 'a25'
    
     # Itero por cada una de las imágenes en el ImageFileCollection
    for hdu, fname in imagenes.hdus(imagetyp='object',return_fname=True):
        
        # Saco informacion del header
        wcs = WCS(hdu.header)
        mjd = hdu.header['mjd-obs']
        gain = hdu.header['gain']
        date = hdu.header['date-obs']
        
        # Imagen
        data = hdu.data

        # calculo la media, desviación y numero de pixeles del background
        bkg_mean = []
        bkg_std = []
        bkg_nsky = []
        anullus_masks = anillos.to_pixel(wcs).to_mask(method='center')
        for mask in anullus_masks:
            annulus_data = mask.multiply(data)
            annulus_data_1d = annulus_data[mask.data > 0]
            nsky = annulus_data_1d.shape[0]
            mean, _, std = sigma_clipped_stats(annulus_data_1d, sigma=2)
            bkg_mean.append(mean)
            bkg_std.append(std)
            bkg_nsky.append(nsky)
        bkg_mean = np.array(bkg_mean)
        bkg_std = np.array(bkg_std)
        bkg_nsky = np.array(bkg_nsky)

        # Área de la apertura
        area=aperturas.to_pixel(wcs).area()

        # realizo fotometria y guardo en tablas ecsv en dir_tablas
        area=aperturas.to_pixel(wcs).area()
        phot = aperture_photometry(data, aperturas, wcs=wcs)
        phot['aper_bkg'] = bkg_mean*area
        phot['aper_sum_bkgsub']=phot['aperture_sum']-phot['aper_bkg']
        inst_mag=[]
        for aper_sum in phot['aper_sum_bkgsub']:
            if aper_sum > 0:
                inst_mag.append(zmag - 2.5*np.log10(aper_sum))
            else:
                inst_mag.append(np.nan)
        phot['inst_mag'] = inst_mag
        err = error(phot['aperture_sum'], gain, area, bkg_std, bkg_nsky)
        phot['inst_mag_err'] = 1.0857 *  err/phot['aperture_sum']
        for col in phot.colnames:
            if col != 'celestial_center':
                phot[col].info.format = '%.8g'
        phot.meta['filename'] = fname  # agrego al metadata el nombre de la imagen
        phot.meta['mjd'] = mjd         # y el mjd
        phot.write(dir_tablas+'/'+fname[:-5]+'.ecsv',overwrite=True)


        # Agrego fila a tabla de fotometría diferencial 
        new_row = [fname, mjd, date]
        for row in phot:
            mag=row['inst_mag']
            merr=row['inst_mag_err']
            new_row.append(mag)
            new_row.append(merr)
        for col in fot_tabla.colnames:
            if col != 'FNAME' and col != 'DATE' and col != 'MJD':
                fot_tabla[col].info.format = '%.8g'

        fot_tabla.add_row(new_row)
    
    clear_output(wait=True)
        
    return fot_tabla



def CalcMejores(tabla, frac=0.1, n=10):    
    '''     
    Está función devuelve las n estrellas más constantes de una
    tabla de fotometría, con sus respectivos pesos
    
    Parámetros
    ----------
    tabla: `~pandas.core.frame.DataFrame` 
        Tabla de la fotometria con sus  respectivos errores las columnas de 
        magnitudes deben tener el nombre de la forma 's#' y los errores 'merr#'
    frac: int, opcional
        Fracción de el total de estrellas que se descarta en cada iteración
    n: int, opcional 
        Numero estrellas que se quiere obtener
        
    Devuelve
    --------
    pesos_finales: `~.pandas.core.series.Series` 
        Fuentes mas constantes del campo con sus respectivos pesos,
        relacionados con el nive de confiabilidad
    '''
    
    
    
    # Armo los indexs para las magnitudes y para los errores
    mags = tabla.columns[tabla.columns.str.contains('s')]
    merrs = tabla.columns[tabla.columns.str.contains('merr')]
    
    # Calculo los pesos iniciales referidos al error en la magnitud
    # los pesos son la inversa de los errores al cuadrado
    # normalizo los pesos
    # cambio el index para que coincida con los de las magnitudes
    pesos_err = tabla[merrs].agg(lambda x: 1/np.mean(x)**2)     
    pesos_err = pesos_err/np.sum(pesos_err)                       
    pesos_err = pesos_err.rename(index = lambda s: 's'+s[4:]) 
    
    # Calculo la dispersión de cada fuente usando los pesos_err
    # pesos_err es un pandas.series
    sigmas = pd.Series(index=mags) # pandas.series de las desviaciones

    for magi in mags:
        
        # Magnitud media de cada imagen
        mag_media = (tabla[mags].agg(axis='columns',
                              func = lambda x : x*pesos_err)
                 .sum(axis='columns'))
    
        # Calculo la desviacion estandar de la fotometría diferencial
        sigma = np.std(tabla[magi] - mag_media)
        # La agrego a "sigmas"
        sigmas[magi] = sigma
        
    ''' Ahora hago un proceso iterativo para calcular las fuentes mas
    constantes y sus respectivos pesos, que seran la inversa de las
    varianzas de sus curvas de luz'''
    mags_mejores = mags
    while len(mags_mejores) > n:     
        # Calculo la curva de luz para cada fuente sin usar dicha fuente
        for magi in mags_mejores:
            
            # Armo el pandas.series sin magi para calcularle su dispersión
            mejores_i = mags_mejores.drop(magi)
            
            # Los pesos los calculo con las dispersiones
            pesos = 1/sigmas[mejores_i]**2
            pesos = pesos / np.sum(pesos)
            
            # Magnitud media de cada imagen
            mag_media = (tabla[mejores_i].agg(axis='columns',
                                          func = lambda x : x*pesos)
                         .sum(axis='columns'))
            
            # Calculo las desviaciones estandar de la fotometría diferencia
            sigma = np.std(tabla[magi] - mag_media)
            # Las agrego a el pandas.series "sigmas"
            sigmas[magi] = sigma
            
        # Quito de la lista la peor fracción 
        descarte = - round(frac * len(mags_mejores))
        if (len(mags_mejores) + descarte) < n or descarte == 0:
            descarte = -1
                
        # La próxima iteración será sobre menos fuentes
        mags_mejores = sigmas.sort_values().index[:descarte]
        sigmas = sigmas[mags_mejores]

    pesos_finales = 1/sigmas**2
    pesos_finales = pesos_finales/pesos_finales.sum()
    
    return pesos_finales