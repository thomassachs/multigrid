import numpy as np
import matplotlib.pyplot as plt
import timeit
from mpl_toolkits.mplot3d import Axes3D
from functions import *


""" Startwerte """
dimension_des_gitters = n = 65      # sollte die Form 2^x +1 haben
homogene_randbedingung = 0
grundgebietstart_x = 0
grundgebietstart_y = 0
grundgebiet_breite = 1
iterationswiederholung = 1
relaxationsfaktor = 1
maximale_zyklenzahl = 10000
konvergenzkriterium = 0.01
plot_nach_jedem_zyklus = False
                
""" rechte seite erstellen """
rechte_seite = np.zeros((n,n))
for i in range(1, n-1):
        for j in range(1, n-1):
            
            x = grundgebietstart_x + i * grundgebiet_breite/(n-1)
            y = grundgebietstart_y + j * grundgebiet_breite/(n-1)
            
            # Funktionswert multipliziert mit der Gitterweite^2
            rechte_seite[i,j] = (200*((1-6*x**2)*y**2*(1-y**2) + (1-6*y**2) * x**2*(1-x**2))
                        * (grundgebiet_breite / (n-1))**2) 

""" exakte lsg erstellen """
exakte_lsg = np.zeros((n,n))
  
for i in range(1, n-1):
    for j in range(1, n-1):
        x = grundgebietstart_x + i * grundgebiet_breite/(n-1)
        y = grundgebietstart_y + j * grundgebiet_breite/(n-1)
        
        exakte_lsg[i, j] = 100*(x**2 - x**4)*(y**4-y**2)
    
""" Interpolations-, Restriktions- und Relaxationsmethode festlegen """
def interpolation(defekt, homogene_randbedingung):
    # momentan nur full_weighting_prolongation
    return full_weighting_prolongation(defekt,homogene_randbedingung)

def restriktion(residuum, homogene_randbedingung):
    # Auswahl zwischen
    # full_weighting_restriktion
    # half_weighting_restriktion
    # injektion_restriktion
    return injektion_restriktion(residuum, homogene_randbedingung)
    
def relaxation(feines_gitter,grundgebietstart_x,grundgebietstart_y,grundgebiet_breite, rechte_seite, relaxationsfaktor):
    # Auswahl zwischen
    # gesamtschritt_lexikographisch
    # einzelschritt_lexikographisch
    return gesamtschritt_lexikographisch(feines_gitter,grundgebietstart_x,grundgebietstart_y,grundgebiet_breite, rechte_seite, relaxationsfaktor)











""" kein gitter nur Iterationsmethode """

""" Zeit stoppen """
startzeit = timeit.default_timer()

# feines Gitter erstellen
feines_gitter = gitter_mit_randbedingung_erstellen(dimension_des_gitters,homogene_randbedingung)

for y in range(maximale_zyklenzahl):
    # altes gitter merken
    altes_gitter = feines_gitter
    
    # feines Gitter iterieren
    feines_gitter = relaxation(feines_gitter,grundgebietstart_x,grundgebietstart_y,grundgebiet_breite, rechte_seite, relaxationsfaktor)
    
    if (plot_nach_jedem_zyklus):
        plot(feines_gitter,grundgebietstart_x, grundgebietstart_y, grundgebiet_breite, exakte_lsg, y)
        
    if (np.amax(abs(exakte_lsg - feines_gitter)) < konvergenzkriterium):
        print('Konvergenzkriterium nach ', y, ' Iterationen erreicht')
        break

# plot nach dem wir mit das Konvergenzkriterium geschafft haben oder die maximale Zykluszahl erreicht wurde    
#plot(feines_gitter,grundgebietstart_x, grundgebietstart_y, grundgebiet_breite, exakte_lsg)

endzeit = timeit.default_timer()

print("runtime Keingitter: {} s".format(float(round(endzeit - startzeit, 3))))






""" Zweigitter-Zyklus """

""" Zeit stoppen """
startzeit = timeit.default_timer()

# feines Gitter erstellen
feines_gitter = gitter_mit_randbedingung_erstellen(dimension_des_gitters,homogene_randbedingung)

for y in range(maximale_zyklenzahl):
    # altes gitter merken
    altes_gitter = feines_gitter
    # feines Gitter iterieren
    for x in range(iterationswiederholung):
        feines_gitter = relaxation(feines_gitter,grundgebietstart_x,grundgebietstart_y,grundgebiet_breite, rechte_seite, relaxationsfaktor)

    # residuum berechnen
    residuum = (residuum_berechnen(feines_gitter,grundgebietstart_x,grundgebietstart_y,grundgebiet_breite,rechte_seite))

    # residuum auf das grobe Gitter restringieren
    residuum_restringiert = restriktion(residuum, homogene_randbedingung)

    # defektgleichung mit restringiertem residuum als rechter Seite iterieren
    defekt = np.zeros((residuum_restringiert.shape))
    for x in range(iterationswiederholung):
        defekt = relaxation(defekt,grundgebietstart_x,grundgebietstart_y,grundgebiet_breite, residuum_restringiert, relaxationsfaktor)

    # defekt prolongieren
    prolongierter_defekt = interpolation(defekt,homogene_randbedingung)

    # korrektur des feinen Gitters
    feines_gitter = feines_gitter + prolongierter_defekt
    
    # feines Gitter iterieren
    for x in range(iterationswiederholung):
        feines_gitter = relaxation(feines_gitter,grundgebietstart_x,grundgebietstart_y,grundgebiet_breite, rechte_seite, relaxationsfaktor)
    
    if (plot_nach_jedem_zyklus):
        plot(feines_gitter,grundgebietstart_x, grundgebietstart_y, grundgebiet_breite, exakte_lsg, y)
        
    if (np.amax(abs(exakte_lsg - feines_gitter)) < konvergenzkriterium):
        print('Konvergenzkriterium nach ', y, ' Zyklen erreicht')
        break
    
# plot nach dem wir mit das Konvergenzkriterium geschafft haben oder die maximale Zykluszahl erreicht wurde    
#plot(feines_gitter,grundgebietstart_x, grundgebietstart_y, grundgebiet_breite, exakte_lsg)

endzeit = timeit.default_timer()

print("runtime Zweigitter: {} s".format(float(round(endzeit - startzeit, 3))))





""" Dreigitter V-Zyklus """


""" Zeit stoppen """
startzeit = timeit.default_timer()
# die Zahlen hinter dem Residuum etc. stehen für das Gitter auf dem wir uns befinden
# 0 steht also für das feinste, 1 für das grobere und 2 für das gröbste Gitter

# feines Gitter erstellen
feines_gitter = gitter_mit_randbedingung_erstellen(dimension_des_gitters,homogene_randbedingung)

for y in range(maximale_zyklenzahl):
    
    # altes gitter merken
    altes_gitter = feines_gitter
    # feines Gitter iterieren
    for x in range(iterationswiederholung):
        feines_gitter = relaxation(feines_gitter,grundgebietstart_x,grundgebietstart_y,grundgebiet_breite, rechte_seite, relaxationsfaktor)


    # residuum0 berechnen
    residuum0 = (residuum_berechnen(feines_gitter,grundgebietstart_x,grundgebietstart_y,grundgebiet_breite, rechte_seite))

    # residuum0 auf das grobe Gitter restringieren
    residuum_restringiert1 = restriktion(residuum0, homogene_randbedingung)

    # defektgleichung mit restringiertem residuum1 als rechter Seite iterieren 
    defekt1 = np.zeros((residuum_restringiert1.shape))
    for x in range(iterationswiederholung):
        defekt1 = relaxation(defekt1,grundgebietstart_x,grundgebietstart_y,grundgebiet_breite, residuum_restringiert1, relaxationsfaktor)

    
    # residuum1 berechnen
    residuum1 = (residuum_berechnen(defekt1,grundgebietstart_x,grundgebietstart_y,grundgebiet_breite, residuum_restringiert1))

    # residuum1 auf das grobe Gitter restringieren
    residuum_restringiert2 = restriktion(residuum1, homogene_randbedingung)

    # defektgleichung mit restringiertem residuum2 als rechter Seite iterieren 
    defekt2 = np.zeros((residuum_restringiert2.shape))
    for x in range(iterationswiederholung):
        defekt2 = relaxation(defekt2,grundgebietstart_x,grundgebietstart_y,grundgebiet_breite, residuum_restringiert2, relaxationsfaktor)

    # defekt2 prolongieren
    prolongierter_defekt1 = interpolation(defekt2,homogene_randbedingung)

    # korrektur des zweiten Gitters
    defekt1 = defekt1 + prolongierter_defekt1

    # grobes Gitter iterieren
    for x in range(iterationswiederholung):
        defekt1 = relaxation(defekt1,grundgebietstart_x,grundgebietstart_y,grundgebiet_breite, residuum_restringiert1, relaxationsfaktor)

    # defekt1 prolongieren
    prolongierter_defekt0 = interpolation(defekt1,homogene_randbedingung)

    # korrektur des feinen Gitters
    feines_gitter = feines_gitter + prolongierter_defekt0

    # feines Gitter iterieren
    for x in range(iterationswiederholung):
        feines_gitter = relaxation(feines_gitter,grundgebietstart_x,grundgebietstart_y,grundgebiet_breite, rechte_seite, relaxationsfaktor)
    
    if (plot_nach_jedem_zyklus):
        plot(feines_gitter,grundgebietstart_x, grundgebietstart_y, grundgebiet_breite, exakte_lsg)
        
    if (np.amax(abs(exakte_lsg - feines_gitter)) < konvergenzkriterium):
        print('Konvergenzkriterium nach ', y, ' Zyklen erreicht')
        break

# plot nach dem wir mit das Konvergenzkriterium geschafft haben oder die maximale Zykluszahl erreicht wurde
#plot(feines_gitter,grundgebietstart_x, grundgebietstart_y, grundgebiet_breite, exakte_lsg)


endzeit = timeit.default_timer()

print("runtime Dreigitterzyklus: {} s".format(float(round(endzeit - startzeit, 3))))




