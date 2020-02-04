#!/usr/bin/python
# -*- coding: utf-8 -*-


"""
Created on Wed Oct 24 13:21:04 2018

@author: eleonoraalei
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
import pyvo
from guizero import CheckBox, App, PushButton,  info, error, \
    Text, Window, Picture, ListBox, TextBox
import matplotlib.pyplot as plt
import time
import os
from matplotlib.lines import Line2D
import fnmatch
from PIL import Image

service = pyvo.dal.TAPService("http://archives.ia2.inaf.it/vo/tap/projects")
response = service.run_sync("SELECT * from exomercat.exomercat",timeout=None)
table=response.to_table()
catalog=table.to_pandas()
if len(catalog)>0:
        for col in catalog.columns:
          try:
            catalog[col]=catalog.apply(lambda row: str(row[col]).lstrip('b').strip("'").replace('"',''),axis=1)
          except:
              pass

for c in catalog.columns:
    try:
        catalog[c]=catalog[c].astype(float)
    except:
        pass
catalog = catalog.replace('nan', np.nan)



# %%% UTILITIES

def is_float(value):
    try:
        float(value)
        return True
    except:
        return False


def write_update():


    note = plt.annotate(
        'LAST UPDATE ' + time.strftime('%d/%m/%Y'),
        xy=(0.7, 0.),
        xytext=(0, 8),
        xycoords=('axes fraction', 'figure fraction'),
        textcoords='offset points',
        size=6,
        ha='left',
        va='bottom',
        )

    return note


def enablelist():
    if listcheck.value == 1:
        listbox.enabled = False
    else:

        listbox.enabled = True


# %%

def reduce_catalog(
    catalog,
    massmin,
    massmax,
    masscheck,
    msinicheck,
    koicheck,
    radmin,
    radmax,
    permin,
    permax,
    amin,
    amax,
    emin,
    emax,
    imin,
    imax,
    listbox,
    folder,
    ):

    if masscheck == 1 and msinicheck == 0:
        if is_float(massmin):
            catalog = catalog[catalog['mass'] >= float(massmin)]
        if is_float(massmax):
            catalog = catalog[catalog['mass'] <= float(massmax)]
    elif masscheck == 0 and msinicheck == 1:
        if is_float(massmin):
            catalog = catalog[catalog['msini'] >= float(massmin)]
        if is_float(massmax):
            catalog = catalog[catalog['msini'] <= float(massmax)]
    else:
        if is_float(massmin):
            catalog = catalog[catalog['bestmass'] >= float(massmin)]
        if is_float(massmax):
            catalog = catalog[catalog['bestmass'] <= float(massmax)]
    if is_float(radmin):
        catalog = catalog[catalog['r'] >= float(radmin)]
    if is_float(radmax):
        catalog = catalog[catalog['r'] < float(radmax)]
    if is_float(permin):
        catalog = catalog[catalog['p'] >= float(permin)]
    if is_float(permax):
        catalog = catalog[catalog['p'] < float(permax)]

    if is_float(amin):
        catalog = catalog[catalog['a'] >= float(amin)]
    if is_float(amax):
        catalog = catalog[catalog['a'] < float(amax)]
    if is_float(emin):
        catalog = catalog[catalog['e'] >= float(emin)]
    if is_float(emax):
        catalog = catalog[catalog['e'] < float(emax)]
    if is_float(imin):
        catalog = catalog[catalog['i'] >= float(imin)]
    if is_float(imax):
        catalog = catalog[catalog['i'] < float(imax)]
    if koicheck == 1:
        catalog = catalog[catalog['status'] == 'CONFIRMED']

    catalog = catalog[catalog['discovery_method'].isin(listbox)]
    if not os.path.exists(folder):
        os.makedirs(folder)
    catalog.to_csv(folder + '/' + 'catred.csv')
    return catalog


def plot(catalog, plots, histplots):


    catalog['Colors'] = catalog['discovery_method'].replace({
        'Transit': 'green',
        'Radial Velocity': 'red',
        'Astrometry': 'yellow',
        'Imaging': 'blue',
        'Microlensing': 'orange',
        'Pulsar Timing': 'brown',
        'TTV': 'gray',
        'Other': 'black',
        })
    catalog['Markers'] = catalog.mass_prov.replace({'Mass': 'o',
            'Msini': 'd'})

    for p in plots:

        (fig, ax) = plt.subplots()
        legend_elements = []

        if p[11]:
            plt.errorbar(
                catalog[p[0]],
                catalog[p[2]],
                xerr=[catalog[p[4]], catalog[p[5]]],
                yerr=[catalog[p[6]], catalog[p[7]]],
                ecolor='lightgrey',
                elinewidth=1,
                fmt='+',
                zorder=0,
                label='',
                )
        for (discmeth, group) in catalog.groupby(['discovery_method']):
            col = group.Colors.unique()[0]

            legend_elements.append(Line2D(
                [0],
                [0],
                marker='o',
                color='w',
                label=discmeth,
                markerfacecolor=col,
                markersize=7,
                markeredgecolor='black',
                ))
            for (mass_prov, subgroup) in group.groupby(['Markers']):
                ax.scatter(
                    subgroup[p[0]],
                    subgroup[p[2]],
                    color=col,
                    marker=mass_prov,
                    label=discmeth,
                    edgecolor='black',
                    )

        if p[2] == 'e':
            ax.set(ylim=[0, 1])

        ax.set_xlabel(p[1])
        ax.set_ylabel(p[3])
        if p[8]:
            plt.xscale('log')
        if p[9]:
            plt.yscale('log')
        ax.tick_params(
            axis='both',
            direction='in',
            which='both',
            bottom=True,
            top=True,
            left=True,
            right=True,
            )

        ax.legend(
            handles=legend_elements,
            loc='upper center',
            ncol=4,
            fancybox=True,
            shadow=True,
            bbox_to_anchor=(0.5, 1.20),
            )
        note = write_update()
        plt.tight_layout(pad=2)
        plt.savefig(folder.value + '/' + p[10],
                    bbox_extra_artists=note, dpi=300)
        plt.close()

#        HISTOGRAMS

    for h in histplots:

        (fig, ax) = plt.subplots()
        legend_elements = []
        x = []
        col = []
        for (discmeth, group) in catalog.groupby(['discovery_method']):

            col.append(group.Colors.unique()[0])
            x.append(group[h[0]].dropna())
            legend_elements.append(Line2D(
                [0],
                [0],
                marker='o',
                color='w',
                label=discmeth,
                markerfacecolor=group.Colors.unique()[0],
                markersize=7,
                markeredgecolor='black',
                ))
        binwidth = (max(catalog[h[0]].dropna())
                    - min(catalog[h[0]].dropna())) / 50
        bins = np.arange(min(catalog[h[0]].dropna()),
                         max(catalog[h[0]].dropna()), binwidth)

        ax.legend(
            handles=legend_elements,
            loc='upper center',
            ncol=4,
            fancybox=True,
            shadow=True,
            bbox_to_anchor=(0.5, 1.20),
            )
        ax.set_xlabel(h[1])
        if h[2] == 1:
            if bins[0] == 0:
                plt.hist(x, bins, color=col, edgecolor='black',
                         stacked=True)
            else:
                logbins = np.logspace(np.log10(bins[0]),
                        np.log10(bins[-1]), len(bins))
                plt.xscale('log')
                plt.hist(x, logbins, color=col, edgecolor='black',
                         stacked=True)
        else:

            plt.hist(x, bins, color=col, edgecolor='black',
                     stacked=True)
        if h[3] == 1:
            plt.yscale('log')
        plt.tight_layout(pad=2)
        note = write_update()
        plt.savefig(folder.value + '/' + h[4], bbox_extra_artists=note,
                    dpi=300)
        plt.close()


def check():
    massmax.text_color = 'black'
    massmin.text_color = 'black'
    radmax.text_color = 'black'
    radmin.text_color = 'black'
    permax.text_color = 'black'
    permin.text_color = 'black'
    emax.text_color = 'black'
    emin.text_color = 'black'
    amax.text_color = 'black'
    amin.text_color = 'black'
    imax.text_color = 'black'
    imin.text_color = 'black'
    errortext = ''
    if massmax.value == '':
        massmax.value = 'Any'
    if massmin.value == '':
        massmin.value = 'Any'
    if radmax.value == '':
        radmax.value = 'Any'
    if radmin.value == '':
        radmin.value = 'Any'
    if permax.value == '':
        permax.value = 'Any'
    if permin.value == '':
        permin.value = 'Any'
    if emax.value == '':
        emax.value = 'Any'
    if emin.value == '':
        emin.value = 'Any'
    if amax.value == '':
        amax.value = 'Any'
    if amin.value == '':
        amin.value = 'Any'
    if imax.value == '':
        imax.value = 'Any'
    if imin.value == '':
        imin.value = 'Any'
    if is_float(massmin.value) and is_float(massmax.value):
        if float(massmin.value) > float(massmax.value):
            massmin.text_color = 'red'
            massmax.text_color = 'red'
            errortext = errortext \
                + 'MASS Measurement error: Minimum is greater than maximum\n'
    else:

        if not is_float(massmin.value) and not massmin.value == 'Any':
            massmin.text_color = 'red'
            errortext = errortext \
                + '''
MASS Type error: Minimum must be a float
'''

        if not is_float(massmax.value) and not massmax.value == 'Any':
            massmax.text_color = 'red'
            errortext = errortext \
                + '''
MASS Type error: Maxmimum must be a float
'''

    if is_float(massmin.value) and float(massmin.value) < 0:
        massmin.text_color = 'red'
        errortext = errortext \
            + '''
MASS Sign error: Minimum must be a positive value
'''

    if is_float(massmax.value) and float(massmax.value) < 0:
        massmax.text_color = 'red'
        errortext = errortext \
            + '''
MASS Sign error: Maximum must be a positive value
'''

    if is_float(radmin.value) and is_float(radmax.value):
        if float(radmin.value) > float(radmax.value):
            radmin.text_color = 'red'
            radmax.text_color = 'red'
            errortext = errortext \
                + '''
RAD Measurement error: Minimum is greater than maximum
'''
    else:
        if not is_float(radmin.value) and not radmin.value == 'Any':
            radmin.text_color = 'red'
            errortext = errortext \
                + '''
RAD Type error: Minimum must be a float
'''

        if not is_float(radmax.value) and not radmax.value == 'Any':
            radmax.text_color = 'red'
            errortext = errortext \
                + '''
RAD Type error: Maxmimum must be a float
'''

    if is_float(radmin.value) and float(radmin.value) < 0:
        radmin.text_color = 'red'
        errortext = errortext \
            + '''
RAD Sign error: Minimum must be a positive value
'''

    if is_float(radmax.value) and float(radmax.value) < 0:
        radmax.text_color = 'red'
        errortext = errortext \
            + '''
RAD Sign error: Maximum must be a positive value
'''

    if is_float(permin.value) and is_float(permax.value):
        if float(permin.value) > float(permax.value):
            permin.text_color = 'red'
            permax.text_color = 'red'
            errortext = errortext \
                + '''
PER Measurement error: Minimum is greater than maximum
'''
    else:

        if not is_float(permin.value) and not permin.value == 'Any':
            permin.text_color = 'red'
            errortext = errortext \
                + '''
PER Type error: Minimum must be a float
'''

        if not is_float(permax.value) and not permax.value == 'Any':
            permax.text_color = 'red'
            errortext = errortext \
                + '''
PER Type error: Maxmimum must be a float
'''
    if is_float(permin.value) and float(permin.value) < 0:
        permin.text_color = 'red'
        errortext = errortext \
            + '''
PER Sign error: Minimum must be a positive value
'''

    if is_float(permax.value) and float(permax.value) < 0:
        permax.text_color = 'red'
        errortext = errortext \
            + '''
PER Sign error: Maximum must be a positive value
'''

    if is_float(amin.value) and is_float(amax.value):
        if float(amin.value) > float(amax.value):
            amin.text_color = 'red'
            amax.text_color = 'red'
            errortext = errortext \
                + '''
A Measurement error: Minimum is greater than maximum
'''
    else:
        if not is_float(amin.value) and not amin.value == 'Any':
            amin.text_color = 'red'
            errortext = errortext \
                + '''
A Type error: Minimum must be a float
'''

        if not is_float(amax.value) and not amax.value == 'Any':
            amax.text_color = 'red'
            errortext = errortext \
                + '''
A Type error: Maxmimum must be a float
'''

    if is_float(amin.value) and float(amin.value) < 0:
        amin.text_color = 'red'
        errortext = errortext \
            + '''
A Sign error: Minimum must be a positive value
'''

    if is_float(amax.value) and float(amax.value) < 0:
        amax.text_color = 'red'
        errortext = errortext \
            + '''
A Sign error: Maximum must be a positive value
'''

    if is_float(emin.value) and is_float(emax.value):
        if is_float(emin.value) and float(emin.value) > 1:
            emin.text_color = 'red'
            errortext = errortext \
                + '''
ECC Value error: Eccentricity must be no greater than 1
'''

        if is_float(emax.value) and float(emax.value) > 1:
            emax.text_color = 'red'
            errortext = errortext \
                + '''
ECC Value error: Eccentricity must be no greater than 1
'''

        if float(emin.value) > float(emax.value):
            emin.text_color = 'red'
            emax.text_color = 'red'
            errortext = errortext \
                + '''
ECC Measurement error: Minimum is greater than maximum
'''
    else:

        if is_float(emin.value) and float(emin.value) > 1:
            emin.text_color = 'red'
            errortext = errortext \
                + '''
ECC Value error: Eccentricity must be no greater than 1
'''

        if is_float(emax.value) and float(emax.value) > 1:
            emax.text_color = 'red'
            errortext = errortext \
                + '''
ECC Value error: Eccentricity must be no greater than 1
'''

        if not is_float(emin.value) and not emin.value == 'Any':
            emin.text_color = 'red'
            errortext = errortext \
                + '''
ECC Type error: Minimum must be a float
'''

        if not is_float(emax.value) and not emax.value == 'Any':
            emax.text_color = 'red'
            errortext = errortext \
                + '''
ECC Type error: Maxmimum must be a float
'''

    if is_float(emin.value) and float(emin.value) < 0:
        emin.text_color = 'red'
        errortext = errortext \
            + '''
ECC Sign error: Minimum must be a positive value
'''

    if is_float(emax.value) and float(emax.value) < 0:
        emax.text_color = 'red'
        errortext = errortext \
            + '''
ECC Sign error: Maximum must be a positive value
'''

    if is_float(emin.value) and float(emin.value) < 0:
        emin.text_color = 'red'
        errortext = errortext \
            + '''
ECC Sign error: Minimum must be a positive value
'''

    if is_float(emax.value) and float(emax.value) < 0:
        emax.text_color = 'red'
        errortext = errortext \
            + '''
ECC Sign error: Maximum must be a positive value
'''

    if is_float(imin.value) and is_float(imax.value):
        if float(imin.value) > float(imax.value):
            imin.text_color = 'red'
            imax.text_color = 'red'
            errortext = errortext \
                + '''
I Measurement error: Minimum is greater than maximum
'''
    else:

        if not is_float(imin.value) and not imin.value == 'Any':
            imin.text_color = 'red'
            errortext = errortext \
                + '''
I Type error: Minimum must be a float
'''

        if not is_float(imax.value) and not imax.value == 'Any':
            imax.text_color = 'red'
            errortext = errortext \
                + '''
I Type error: Maxmimum must be a float
'''

    if is_float(imin.value) and float(imin.value) < 0:
        imin.text_color = 'red'
        errortext = errortext \
            + '''
I Sign error: Minimum must be a positive value
'''

    if is_float(imax.value) and float(imax.value) < 0:
        imax.text_color = 'red'
        errortext = errortext \
            + '''
I Sign error: Maximum must be a positive value
'''

    if len(errortext) > 0:
        error('Errors', errortext)
    else:
        testo = 'You selected:\nMASS: From ' + str(massmin.value) \
            + ' to ' + str(massmax.value) + '\nMASSCHECK: ' \
            + str(masscheck.value) + '\nMSINICHECK ' \
            + str(msinicheck.value) + '\nRADIUS: From ' \
            + str(radmin.value) + ' to ' + str(radmax.value) \
            + '\nPERIOD: From ' + str(permin.value) + ' to ' \
            + str(permax.value) + '\nSEMI-MAJOR AXIS: From ' \
            + str(amin.value) + ' to ' + str(amax.value) \
            + '\nECCENTRICITY: From ' + str(emin.value) + ' to ' \
            + str(emax.value) + '\nINCLINATION: From ' \
            + str(imin.value) + ' to ' + str(imax.value) \
            + '\nDiscovery methods: ' + str(listbox.value)
        info('Recap', testo)

    return errortext


def stdplot():
    error = check()
    if len(error) > 0:
        pass
    else:
        catred = reduce_catalog(
            catalog,
            massmin.value,
            massmax.value,
            masscheck.value,
            msinicheck.value,
            koi.value,
            radmin.value,
            radmax.value,
            permin.value,
            permax.value,
            amin.value,
            amax.value,
            emin.value,
            emax.value,
            imin.value,
            imax.value,
            listbox.value,
            folder.value,
            )

        plots = [[
            'a',
            'Semi-major axis (AU)',
            'bestmass',
            r'Best Mass (M$_J$)',
            'a_max',
            'a_min',
            'bestmass_max',
            'bestmass_min',
            1,
            1,
            'BMASSa.png',
            1,
            ], [
            'bestmass',
            r'Best Mass (M$_J$)',
            'r',
            'Radius (R$_J$)',
            'bestmass_max',
            'bestmass_min',
            'r_max',
            'r_min',
            1,
            1,
            'BMASSRAD.png',
            1,
            ], [
            'p',
            'Period (days)',
            'bestmass',
            r'Best Mass (M$_J$)',
            'p_max',
            'p_min',
            'bestmass_max',
            'bestmass_min',
            1,
            1,
            'BMASSPER.png',
            1,
            ], [
            'bestmass',
            r'Best Mass (M$_J$)',
            'e',
            r'Eccentricity',
            'bestmass_max',
            'bestmass_min',
            'e_max',
            'e_min',
            1,
            0,
            'BMASSe.png',
            1,
            ], [
            'a',
            'Semi-major axis (AU)',
            'e',
            r'Eccentricity',
            'a_max',
            'a_min',
            'e_max',
            'e_min',
            1,
            0,
            'ea.png',
            1,
            ]]

        histplots = [
            ['a', 'Semi-major axis (AU)', 1, 1, 'hista.png'],
            ['mass', r'Mass (M$_J$)', 1, 1, 'histm.png'],
            ['msini', r'Minimum Mass (M$_J$)', 1, 1, 'histmsini.png'],
            ['r', 'Radius (R$_J$)', 1, 1, 'histr.png'],
            ['p', 'Period (days)', 1, 1, 'histp.png'],
            ['e', r'Eccentricity', 0, 1, 'histe.png'],
            ['i', 'Inclination (deg)', 0, 1, 'histi.png'],
            ]

        plot(catred, plots, histplots)
        show_plots()


def show_plots():

    window = Window(app, title='Plots window', height=1000, width=1300,
                    layout='grid')
    i = 0
    j = 0
    for entry in os.listdir(folder.value):

        if fnmatch.fnmatch(entry, '[!h]*.png'):

            picture = Picture(window, image=folder.value + '/' + entry,
                              grid=[i, j],width=300, height=200)
            i = i + 1
            if i % 4 == 0:
                j = j + 1
                i = 0
    for entry in os.listdir(folder.value):

        if fnmatch.fnmatch(entry, 'h*.png'):

            picture = Picture(window, image=folder.value + '/' + entry,
                              grid=[i, j],width=300, height=200)
            i = i + 1
            if i % 4 == 0:
                j = j + 1
                i = 0


def confirm_settings(
    catred,
    mr,
    mrx,
    mry,
    mp,
    mpx,
    mpy,
    msma,
    msmax,
    msmay,
    me,
    mex,
    mey,
    esma,
    esmax,
    esmay,
    hista,
    histax,
    histay,
    histm,
    histmx,
    histmy,
    histmsini,
    histmsinix,
    histmsiniy,
    histr,
    histrx,
    histry,
    histp,
    histpx,
    histpy,
    histi,
    histiy,
    histe,
    histey,
    mrerr,
    msmaerr,
    mperr,
    esmaerr,
    meerr,
    ):

    plots = []
    if msma.value == 1:
        plots.append([
            'a',
            'Semi-major axis (AU)',
            'bestmass',
            r'Best Mass (M$_J$)',
            'a_max',
            'a_min',
            'bestmass_max',
            'bestmass_min',
            msmax.value,
            msmay.value,
            'BMASSa.png',
            msmaerr.value,
            ])
    if mr.value == 1:
        plots.append([
            'r',
            r'Radius (R$_J$))',
            'bestmass',
            r'Best Mass (M$_J$)',
            'r_max',
            'r_min',
            'bestmass_max',
            'bestmass_min',
            mrx.value,
            mry.value,
            'BMASSRAD.png',
            mrerr.value,
            ])
    if mp.value == 1:
        plots.append([
            'bestmass',
            r'Best Mass (M$_J$)',
            'p',
            'Period (days)',
            'bestmass_max',
            'bestmass_min',
            'p_max',
            'p_min',
            mpx.value,
            mpy.value,
            'BMASSPER.png',
            mperr.value,
            ])
    if me.value == 1:
        plots.append([
            'bestmass',
            r'Best Mass  (M$_J$)',
            'e',
            r'Eccentricity',
            'bestmass_max',
            'bestmass_min',
            'e_max',
            'e_min',
            mex.value,
            mey.value,
            'BMASSe.png',
            meerr.value,
            ])
    if esma.value == 1:
        plots.append([
            'a',
            'Semi-major axis (AU)',
            'e',
            'Eccentricity',
            'a_max',
            'a_min',
            'e_max',
            'e_min',
            esmax.value,
            esmay.value,
            'ea.png',
            esmaerr.value,
            ])

    histplots = []
    if hista.value == 1:
        histplots.append(['a', 'Semi-major axis (AU)', histax.value,
                         histay.value, 'hista.png'])
    if histm.value == 1:
        histplots.append(['mass', r'Mass (M$_J$)', histmx.value,
                         histmy.value, 'histm.png'])
    if histmsini.value == 1:
        histplots.append(['msini', r'Minimum Mass (M$_J$)', histmsinix.value,
                         histmsiniy.value, 'histmsini.png'])
    if histr.value == 1:
        histplots.append(['r', 'Radius (R$_J$)', histrx.value,
                         histry.value, 'histr.png'])
    if histp.value == 1:
        histplots.append(['p', 'Period (days)', histpx.value,
                         histpy.value, 'histp.png'])
    if histe.value == 1:
        histplots.append(['e', r'Eccentricity', 0, histey.value,
                         'histe.png'])
    if histi.value == 1:
        histplots.append(['i', 'Inclination (deg)', 0, histiy.value,
                         'histi.png'])

    plot(catred, plots, histplots)

    show_plots()


def advchoice(catred):
    window2 = Window(app, title='Advanced settings', height=300,
                     width=600, layout='grid')

    mr = CheckBox(window2, text='Radius (x) vs Mass (y)', grid=[0, 0],
                  align='left')
    mr.value = 1
    mrx = CheckBox(window2, text='logscale x', grid=[1, 0], align='left'
                   )
    mrx.value = 1
    mry = CheckBox(window2, text='logscale y', grid=[2, 0], align='left'
                   )
    mry.value = 1
    mrerr = CheckBox(window2, text='errorbar', grid=[3, 0], align='left'
                     )
    mrerr.value = 1
    mp = CheckBox(window2, text='Mass (x) vs Period (y)', grid=[0, 1],
                  align='left')
    mp.value = 1
    mpx = CheckBox(window2, text='logscale x', grid=[1, 1], align='left'
                   )
    mpx.value = 1
    mpy = CheckBox(window2, text='logscale y', grid=[2, 1], align='left'
                   )
    mpy.value = 1
    mperr = CheckBox(window2, text='errorbar', grid=[3, 1], align='left'
                     )
    mperr.value = 1
    msma = CheckBox(window2, text='Semi-major axis (x) vs Mass (y)',
                    grid=[0, 2], align='left')
    msma.value = 1
    msmax = CheckBox(window2, text='logscale x', grid=[1, 2],
                     align='left')
    msmax.value = 1
    msmay = CheckBox(window2, text='logscale y', grid=[2, 2],
                     align='left')
    msmay.value = 0
    msmaerr = CheckBox(window2, text='errorbar', grid=[3, 2],
                       align='left')
    msmaerr.value = 1
    me = CheckBox(window2, text='Mass (x) vs Eccentricity (y)',
                  grid=[0, 3], align='left')
    me.value = 1
    mex = CheckBox(window2, text='logscale x', grid=[1, 3], align='left'
                   )
    mex.value = 1
    mey = CheckBox(window2, text='logscale y', grid=[2, 3], align='left'
                   )
    mey.value = 0
    meerr = CheckBox(window2, text='errorbar', grid=[3, 3], align='left'
                     )
    meerr.value = 1
    esma = CheckBox(window2,
                    text='Eccentricity (x) vs Semi-major axis (y)',
                    grid=[0, 4], align='left')
    esma.value = 1
    esmax = CheckBox(window2, text='logscale x', grid=[1, 4],
                     align='left')
    esmax.value = 1
    esmay = CheckBox(window2, text='logscale y', grid=[2, 4],
                     align='left')
    esmay.value = 0
    esmaerr = CheckBox(window2, text='errorbar', grid=[3, 4],
                       align='left')
    esmaerr.value = 1
    hista = CheckBox(window2, text='Histogram Semi-major Axis',
                     grid=[0, 5], align='left')
    hista.value = 1
    histax = CheckBox(window2, text='logscale x', grid=[1, 5],
                      align='left')
    histax.value = 1
    histay = CheckBox(window2, text='logscale y', grid=[2, 5],
                      align='left')
    histay.value = 1
    histm = CheckBox(window2, text='Histogram Mass', grid=[0, 6],
                     align='left')
    histm.value = 1
    histmx = CheckBox(window2, text='logscale x', grid=[1, 6],
                      align='left')
    histmx.value = 1
    histmy = CheckBox(window2, text='logscale y', grid=[2, 6],
                      align='left')
    histmy.value = 1
    histmsini = CheckBox(window2, text='Histogram Msini', grid=[0, 7],
                     align='left')
    histmsini.value = 1
    histmsinix = CheckBox(window2, text='logscale x', grid=[1, 7],
                      align='left')
    histmsinix.value = 1
    histmsiniy = CheckBox(window2, text='logscale y', grid=[2, 7],
                      align='left')
    histmsiniy.value = 1
    histr = CheckBox(window2, text='Histogram Radius', grid=[0, 8],
                     align='left')
    histr.value = 1
    histrx = CheckBox(window2, text='logscale x', grid=[1, 8],
                      align='left')
    histrx.value = 1
    histry = CheckBox(window2, text='logscale y', grid=[2, 8],
                      align='left')
    histry.value = 1
    histp = CheckBox(window2, text='Histogram Period', grid=[0, 9],
                     align='left')
    histp.value = 1
    histpx = CheckBox(window2, text='logscale x', grid=[1, 9],
                      align='left')
    histpx.value = 1
    histpy = CheckBox(window2, text='logscale y', grid=[2, 9],
                      align='left')
    histpy.value = 1
    histe = CheckBox(window2, text='Histogram Eccentricity', grid=[0,
                     10], align='left')
    histe.value = 1
    histey = CheckBox(window2, text='logscale y', grid=[2, 10],
                      align='left')
    histey.value = 1
    histi = CheckBox(window2, text='Histogram Period', grid=[0, 11],
                     align='left')
    histi.value = 1
    histiy = CheckBox(window2, text='logscale y', grid=[2, 11],
                      align='left')
    histiy.value = 1
    button = PushButton(window2, text='Ok', command=confirm_settings,
                        args=[
        catred,
        mr,
        mrx,
        mry,
        mp,
        mpx,
        mpy,
        msma,
        msmax,
        msmay,
        me,
        mex,
        mey,
        esma,
        esmax,
        esmay,
        hista,
        histax,
        histay,
        histm,
        histmx,
        histmy,
        histmsini,
        histmsinix,
        histmsiniy,
        histr,
        histrx,
        histry,
        histp,
        histpx,
        histpy,
        histi,
        histiy,
        histe,
        histey,
        mrerr,
        msmaerr,
        mperr,
        esmaerr,
        meerr,
        ], grid=[4, 13, 2, 1])


def advplot():
    error = check()
    if len(error) > 0:
        pass
    else:
        catred = reduce_catalog(
            catalog,
            massmin.value,
            massmax.value,
            masscheck.value,
            msinicheck.value,
            koi.value,
            radmin.value,
            radmax.value,
            permin.value,
            permax.value,
            amin.value,
            amax.value,
            emin.value,
            emax.value,
            imin.value,
            imax.value,
            listbox.value,
            folder.value,
            )
        advchoice(catred)



# %% MAIN

app = App(title='Main window', layout='grid', height=400, width=600)

text = Text(app, text='MINIMUM', grid=[2, 0])
text.height = 2
text = Text(app, text='MAXIMUM', grid=[3, 0])
text.height = 2
text = Text(app, text='Parameter', grid=[1, 0])
text.width = 15
text = Text(app, text='Unit', grid=[4, 0], align='left')
text.width = 10
text = Text(app, text='Mass ', grid=[1, 2], align='right')
unit = Text(app, text='M_J', grid=[4, 2], align='left')

massmin = TextBox(app, text='', grid=[2, 2])
massmax = TextBox(app, text='', grid=[3, 2])

text = Text(app, text='Radius', grid=[1, 3], align='right')
text.height = 2
unit = Text(app, text='R_J', grid=[4, 3], align='left')
unit.height = 2
radmin = TextBox(app, text='', grid=[2, 3])
radmax = TextBox(app, text='', grid=[3, 3])

text = Text(app, text='Period', grid=[1, 4], align='right')
text.height = 2
unit = Text(app, text='days', grid=[4, 4], align='left')
unit.height = 2
permin = TextBox(app, text='', grid=[2, 4])
permax = TextBox(app, text='', grid=[3, 4])

text = Text(app, text='Semi-major axis', grid=[1, 5], align='right')
text.height = 2
unit = Text(app, text='AU', grid=[4, 5], align='left')
unit.height = 2
amin = TextBox(app, text='', grid=[2, 5])
amax = TextBox(app, text='', grid=[3, 5])

text = Text(app, text=r'Eccentricity', grid=[1, 6], align='right')
text.height = 2
unit = Text(app, text='', grid=[4, 6])
unit.height = 2
emin = TextBox(app, text='', grid=[2, 6])
emax = TextBox(app, text='', grid=[3, 6])

text = Text(app, text='Inclination', grid=[1, 7], align='right')
text.height = 2
unit = Text(app, text='degrees', grid=[4, 7], align='left')
unit.height = 2
imin = TextBox(app, text='', grid=[2, 7])
imax = TextBox(app, text='', grid=[3, 7])

massmin.value = 'Any'
massmax.value = 'Any'
radmin.value = 'Any'
radmax.value = 'Any'
permin.value = 'Any'
permax.value = 'Any'
amin.value = 'Any'
amax.value = 'Any'
emin.value = 'Any'
emax.value = 'Any'
imin.value = 'Any'
imax.value = 'Any'

text = Text(app, text='Folder Name', grid=[0, 8, 2, 1], align='right')
text.height = 2
text.width = 15
folder = TextBox(app, grid=[2, 8, 2, 1])


folder.value = time.strftime('%Y%m%d') + '/'
folder.width = 22

koi = CheckBox(app, text='only confirmed', grid=[9, 0], align='left')
koi.height = 2
masscheck = CheckBox(app, text='Msini', grid=[9, 2], align='left')
masscheck.value = 1
masscheck.height = 2
msinicheck = CheckBox(app, text='Mass', grid=[10, 2], align='left')
msinicheck.value = 1
msinicheck.height = 2

text = Text(app, text='Discovery Method', grid=[9, 3])
listbox = ListBox(
    app,
    items=[
        'Radial Velocity',
        'Transit',
        'Astrometry',
        'Imaging',
        'Microlensing',
        'TTV',
        'Pulsar Timing',
        'Other',
        ],
    multiselect=True,
    grid=[7, 4, 4, 4],
   # width=15,
    #height=8,
    )
listbox.value = [
    'Radial Velocity',
    'Transit',
    'Astrometry',
    'Imaging',
    'Microlensing',
    'TTV',
    'Pulsar Timing',
    'Other',
    ]
listbox.enabled = False
listcheck = CheckBox(app, text='All', grid=[10, 3], align='left',
                     command=enablelist)
listcheck.value = 1

button = PushButton(app, text='Plot', command=stdplot, grid=[4, 11, 2,
                    1])
button = PushButton(app, text='Advanced Plot', command=advplot,
                    grid=[1, 11, 2, 1])
app.display()
