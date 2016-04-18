#!/usr/bin/python2.7
# encoding: utf-8

import os
import datetime
import matplotlib.cm as cmap
import numpy as np
from decimal import Decimal
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import A4, landscape
from reportlab.platypus import BaseDocTemplate, Frame, Paragraph, PageTemplate
from reportlab.platypus import Spacer, Table, Image, NextPageTemplate, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch, cm
from reportlab.lib.utils import ImageReader

# Custom error
from pyseidon_dvt.utilities.pyseidon_error import PyseidonError


def write_report(valClass, report_title="validation_report.pdf", debug=False):
    """
    Write automatically and dynamically a report based on the validation processed

    Args:
      valClass (pyseidon.Validation): Validation object

    Kwargs:
      report_title (str): report file name
      debug (bool): debug flag
    """
    # tests
    if not (hasattr(valClass, "Benchmarks") or hasattr(valClass, "HarmonicBenchmarks")):
        raise PyseidonError("-Run validation function first-")

    benchflag = False
    if hasattr(valClass, "Benchmarks"):
        benchflag = True

    harmoflag = False
    if hasattr(valClass, "HarmonicBenchmarks"):
        harmoflag = True

    # date
    now = datetime.date.today()
    # initiate doc
    doc = BaseDocTemplate(report_title,
                            pagesize=A4,
                            title="Validation Report " + now.strftime("- %d %B %Y"),
                            rightMargin=72,leftMargin=72,
                            topMargin=72,bottomMargin=18)

    # default style and frames
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))
    pgtemp = []
    imNb = -1

    #normal frame as for SimpleFlowDocument
    frameP = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id='normal')  # portray frame
    frameL = Frame(doc.leftMargin * 0.75,  doc.bottomMargin, doc.height, doc.width, id='normal')  # landscape frame

    #Two Columns
    frame1 = Frame(doc.leftMargin, doc.bottomMargin, doc.width/2-6, doc.height, id='col1')
    frame2 = Frame(doc.leftMargin+doc.width/2+6, doc.bottomMargin, doc.width/2-6, doc.height, id='col2')

    # initialize "story" list
    story = []

    # Title
    story.append(Spacer(doc.width, A4[1]/3.0))
    story.append(Paragraph("Validation benchmarks and methodologies for FVCOM Hydrodynamic Model results",
                           styles['Title']))
    story.append(Spacer(doc.width, A4[1]/3.0))
    story.append(Paragraph(now.strftime("%B, %Y"), styles['BodyText']))
    pgtemp.append(PageTemplate(id='OneCol', frames=frameP, onPage=footpagenumber))

    # Introduction
    story.append(NextPageTemplate('OneCol'))
    story.append(PageBreak())
    story.append(Paragraph("Introduction", styles['Heading1']))
    story.append(Paragraph("The following report is a summary of validation standards and \
                           methodologies usually performed on the FVCOM model output data in comparison \
                           to observed data. This collection and analysis of a set of \
                           statistics mostly adhere to the benchmarks defined as standards for \
                           hydrodynamic model validation by NOAA [1]. Additional statistics have been \
                           added to provide additional clarity on the skill of the model [2, 3, 4]."
                           , styles['BodyText']))
    story.append(Paragraph("The present validation set is performed the following variables:"
                           , styles['BodyText']))
    if benchflag:
        story.append(Paragraph("Hydrodynamic quantities", styles['Italic']))
        story.append(Paragraph("el: Elevation (deviation from mean sea level, m)"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("cubic_speed: Signed cubic flow speed (m/s)"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("u: normal velocity component (m/s)"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("v: tangent velocity component (m/s)"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("vel: signed flow velocity (m/s)"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("speed: flow speed (m/s)"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("dir: Current direction (between -180 and 180 degrees)"
                               , styles['Bullet'], bulletText='-'))
    if harmoflag:
        story.append(Paragraph("Harmonic coefficients", styles['Italic']))
        story.append(Paragraph("A & A_ci: Amplitude and associated 95% confidence interval"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("g & g_ci: Greenwich phase lag and associated 95% confidence interval"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("theta & theta_ci: Current ellipse orientation angle and associated 95% confidence interval"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("Lsmin & Lsmin_ci: Current ellipse minor axis length and associated 95% confidence interval"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("Lsmaj & Lsmaj_ci: Current ellipse major axis length and associated 95% confidence interval"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("constituent: tidal constituent's name"
                               , styles['Bullet'], bulletText='-'))
    story.append(Spacer(1, 12))
    if benchflag and harmoflag:
        pgtemp.append(PageTemplate(id='OneCol', frames=frameP, onPage=footpagenumber))

    # Statistics
    if benchflag and harmoflag:
        story.append(NextPageTemplate('OneCol'))
        story.append(PageBreak())

    story.append(Paragraph("Statistics", styles['Heading1']))
    if benchflag:
        story.append(Paragraph("Following is a list of the statistics used to evaluate model skill, separated \
                                into two categories: NOAA's Standard Suite of Statistics (SSS), and those added \
                                (Additional)"
                               , styles['BodyText']))
        story.append(Paragraph("SSS", styles['Heading2']))
        story.append(Paragraph("RMSE: Root Mean Squared Error"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("NRMSE: Normalized Root Mean Squared Error (in %)"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("MSE: Mean Square Error"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("NMSE: Normalized Mean Square Error (in %)"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("SD: Standard Deviation of Error"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("CF(X): Central Frequency; percentage of error values that fall within the \
                                range (-X, X)"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("POF(X): Positive Outlier Frequency; percentage of error values that fall \
                                above X",
                               styles['Bullet'], bulletText='-'))
        story.append(Paragraph("NOF(X): Negative Outlier Frequency; percentage of error values that fall \
                                below -X"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("MDPO(X): Maximum Duration of Positive Outliers; longest number of \
                                minutes during which consecutive errors fall above X"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("MDNO(X): Maximum Duration of Negative Outliers; longest number of \
                                minutes during which consecutive errors fall below X"
                               , styles['Bullet'], bulletText='-'))

        story.append(Paragraph("Additional", styles['Heading2']))
        story.append(Paragraph("Willmott Skill: A measure of model adherence to observed data between \
                                0 and 1, with 0 being absolutely no adherence, and 1 being perfect"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("Phase: The phase shift (minutes) of model data that minimizes RMSE \
                                across a timespan of +/-3hr"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("R2 (i.e. coefficient of determination): Measure of the strength of the linear \
                                correlation between the model data and the observed data between 0 \
                                and 1, with 0 being no correlation, and 1 being perfect correlation"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("bias: bias of the model, a measure of over/under-estimation"
                               , styles['Bullet'], bulletText='-'))
        story.append(Paragraph("Pbias: percent bias between the model and the observed data"
                               , styles['Bullet'], bulletText='-'))

    if harmoflag:
        story.append(Paragraph("The statistics reported in the 'Harmonic Analysis' section, \
                                can be defined as the normalised error (in %) between the observed and simulated \
                                harmonic analysis coefficients described in the 'Introduction' section."
                               , styles['BodyText']))

    pgtemp.append(PageTemplate(id='OneCol', frames=frameP, onPage=footpagenumber))

    # Results
    story.append(NextPageTemplate('OneCol'))
    story.append(PageBreak())
    story.append(Paragraph("Results", styles['Heading1']))
    story.append(Paragraph("The following map displays the location(s) as well as the type(s) of the measurement(s)  \
                           used in this validation report.",
                           styles['BodyText']))
    story.append(Spacer(1, 12))
    # Map: measurement's locations
    imNb += 1
    savename = 'tmp_'+str(imNb)+'_plot.png'
    lonmax = -1.0 * np.inf
    lonmin = np.inf
    latmax = -1.0 * np.inf
    latmin = np.inf
    for ii, coor in enumerate(valClass._coordinates):
        lon = coor[0]
        lat = coor[1]
        if lon > lonmax: lonmax = lon
        if lat > latmax: latmax = lat
        if lon < lonmin: lonmin = lon
        if lat < latmin: latmin = lat
    #  redefine colorbar min/max
    margin = 0.01
    indices = np.where(np.logical_and(
                       np.logical_and(valClass._simulated.Grid.lon[:] < lonmax + margin,
                                      valClass._simulated.Grid.lon[:] > lonmin - margin),
                       np.logical_and(valClass._simulated.Grid.lat[:] < latmax + margin,
                                      valClass._simulated.Grid.lat[:] > latmin - margin)))[0]
    cmax = valClass._simulated.Grid.h[indices].max()
    cmin = valClass._simulated.Grid.h[indices].min()
    valClass._simulated.Plots.colormap_var(valClass._simulated.Grid.h,
                                           title='Bathymetric Map & Measurement location(s)',
                                           cmax=cmax, cmin=cmin, isoline='var', mesh=False)
    #  redefine frame
    valClass._simulated.Plots._ax.set_xlim([lonmin - margin, lonmax + margin])
    valClass._simulated.Plots._ax.set_ylim([latmin - margin, latmax + margin])
    color = cmap.rainbow(np.linspace(0, 1, len(valClass._coordinates)))
    for ii, coor in enumerate(valClass._coordinates):
        lon = coor[0]
        lat = coor[1]
        name = coor[2]
        txt = str(ii)
        valClass._simulated.Plots._ax.scatter(lon, lat, label=name,
                                              lw=2, s=50, color=color[ii])
        # valClass._simulated.Plots._ax.annotate(txt, (lon, lat), size=20)
    valClass._simulated.Plots._ax.legend()
    valClass._simulated.Plots._fig.savefig(savename, format='png', bbox_inches='tight')
    valClass._simulated.Plots._fig.clear()
    # image = Image(savename, width=doc.width, height=doc.height / 1.5)
    # story.append(image)
    story.append(get_image(savename, width=16*cm))

    pgtemp.append(PageTemplate(id='OneCol', frames=frameP, onPage=footpagenumber))

    # table's style
    ts = [('ALIGN', (1,1), (-1,-1), 'CENTER'),
          ('LINEABOVE', (0,0), (-1,0), 1, colors.black),
          ('LINEBELOW', (0,0), (-1,0), 1, colors.black),
          ('FONT', (0,0), (-1,0), 'Times-Bold')]

    if benchflag:
        story.append(NextPageTemplate('landscape'))
        story.append(PageBreak())
        story.append(Paragraph("Validation Benchmarks", styles['Heading2']))
        df = valClass.Benchmarks
        values = df.values.tolist()
        # clean up
        for ii, l in enumerate(values):
            for jj, e in enumerate(l):
                if type(e) == float:
                    values[ii][jj] = float(Decimal("%.2f" % e))
        lista = [df.columns[:,].values.astype(str).tolist()] + values
        table = Table(lista, style=ts)
        story.append(table)
        pgtemp.append(PageTemplate(id='landscape', frames=frameL, onPage=make_landscape))

    if harmoflag:
        story.append(NextPageTemplate('OneCol'))
        story.append(PageBreak())
        if type(valClass.HarmonicBenchmarks.elevation) != str:
            story.append(Paragraph("Harmonic Analysis - Elevation", styles['Heading2']))
            df = valClass.HarmonicBenchmarks.elevation
            values = df.values.tolist()
            # clean up
            for ii, l in enumerate(values):
                for jj, e in enumerate(l):
                    if type(e) == float:
                        values[ii][jj] = float(Decimal("%.2f" % e))
            lista = [df.columns[:,].values.astype(str).tolist()] + values
            table = Table(lista, style=ts)
            story.append(table)
        if type(valClass.HarmonicBenchmarks.velocity) != str:
            story.append(Paragraph("Harmonic Analysis - Velocity", styles['Heading2']))
            df = valClass.HarmonicBenchmarks.velocity
            values = df.values.tolist()
            # clean up
            for ii, l in enumerate(values):
                for jj, e in enumerate(l):
                    if type(e) == float:
                        values[ii][jj] = float(Decimal("%.2f" % e))
            lista = [df.columns[:,].values.astype(str).tolist()] + values
            table = Table(lista, style=ts)
            story.append(table)
        pgtemp.append(PageTemplate(id='OneCol', frames=frameP, onPage=footpagenumber))

    # Taylor diagram
    if benchflag:
        story.append(NextPageTemplate('OneCol'))
        story.append(PageBreak())
        story.append(Paragraph("The Taylor diagram above is a concise statistical summary of how well patterns match each \
                                other in terms of their correlation, their root-mean-square difference and the ratio of \
                                their variances.",
                               styles['BodyText']))
        story.append(Spacer(1, 12))
        imNb += 1
        savename = 'tmp_'+str(imNb)+'_plot.png'
        valClass.taylor_diagram(savepath="./", fname=savename)
        # image = Image(savename, width=doc.width , height=doc.height / 2.25)
        # story.append(image)
        story.append(get_image(savename, width=16*cm))

        pgtemp.append(PageTemplate(id='OneCol', frames=frameP, onPage=footpagenumber))

    # References
    story.append(NextPageTemplate('OneCol'))
    story.append(PageBreak())
    story.append(Paragraph("References", styles['Heading1']))
    story.append(Paragraph("K. W. Hess, T. F. Gross, R. A. Schmalz, J. G. Kelley, F. Aikman and E. Wei, \
                            NOS standards for evaluating operational nowcast and forecst hydrodynamic model systems, \
                            National Oceanic and Atmospheric Administration, Silver Srping, Maryland, 2003.",
                           styles['Italic'], bulletText='[1]'))
    story.append(Paragraph("K. Gunn and C. Stock-Williams, On validating numerical hydrodynamic models of complex \
                            tidal flow, International Journal of Marine Energy, Vols. 3-4, no. Special, pp. 82-97, \
                            2013.",
                           styles['Italic'], bulletText='[2]'))
    story.append(Paragraph("N. Georgas and A. F. Blumberg, Establishing Confidence in Marine Forecast Systems: \
                            The Design and Skill Assessment of the New York Harbor Observation and Prediction System, \
                            Version 3 (NYHOPS v3), in 11th International Conference on Estuarine and Coastal Modeling, \
                            Seattle, Washington, United States, 2010.",
                           styles['Italic'], bulletText='[3]'))
    story.append(Paragraph("Y. Liu, P. MacCready, H. M. Barbara, E. P. Dever, M. Kosro and N. S. Banas, \
                            Evaluation of a coastal ocean circulation model for the Columbia River plume in summer \
                            2004, Journal of Geophysical Research, vol. 114, no. C2, p. 1978–2012, 2009.",
                           styles['Italic'], bulletText='[4]'))

    pgtemp.append(PageTemplate(id='OneCol', frames=frameP, onPage=footpagenumber))

    # start the construction of the pdf
    doc.addPageTemplates(pgtemp)
    doc.build(story)
    # remove tmp plots
    for ii in range(imNb+1):
        savename = 'tmp_'+str(ii)+'_plot.png'
        os.remove(savename)


def footpagenumber(canvas, doc):
    """
    Foot note and page number
    """
    canvas.saveState()
    canvas.setFont('Times-Roman', 9)
    canvas.drawString(A4[0]/2.0, 0.75 * inch, "%d" % doc.page)
    canvas.restoreState()

def make_landscape(canvas,doc):
    """
    Foot note and page number plus landscape page formatting
    """
    canvas.saveState()
    canvas.setPageSize(landscape(A4))
    canvas.setFont('Times-Roman', 9)
    canvas.drawString(A4[1]/2.0, 0.75 * inch, "%d" % doc.page)
    canvas.restoreState()

def get_image(path, width=1*cm):
    """
    Read image from file and fit its size to the given width
    """
    img = ImageReader(path)
    iw, ih = img.getSize()
    aspect = ih / float(iw)
    return Image(path, width=width, height=(width * aspect))