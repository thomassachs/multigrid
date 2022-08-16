import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def startgitter_erstellen(n):
    # anfangslsg erstellen
    # default ist ein null Array n x n
    return np.zeros((n,n))
     #return np.random.random((n, n))

def randbedingung_erstellen(n,homogene_randbedingung=0):
    # randbedingung wird hier definiert
    # Matrix mit n x n, wobei die inneren Punkte 0 gesetzt werden und die umliegenden Punkte mit der homogenen Randbedinung
    
    randbedingung = np.zeros((n,n))
    
    # oben
    randbedingung[0, 0:n] = homogene_randbedingung
    
    #links
    randbedingung[0:n, 0] = homogene_randbedingung
    
    #rechts
    randbedingung[0:n, n-1] = homogene_randbedingung
    
    #unten
    randbedingung[n-1, 0:n] = homogene_randbedingung
    
    return randbedingung

def gitter_mit_randbedingung_erstellen(n, homogene_randbedingung=0):
    
    startgitter = startgitter_erstellen(n-2)
    randbedingung = randbedingung_erstellen(n,homogene_randbedingung)
    randbedingung[1:n-1, 1:n-1] = startgitter
    
    return randbedingung

def gesamtschritt_lexikographisch(T,grundgebietstart_x, grundgebietstart_y, grundgebiet_breite, rechte_seite, relaxationsfaktor):

    m, n = T.shape

    _T = np.copy(T)

    # iterate over interior

    for i in range(1, m-1):
        for j in range(1, n-1):
            
            x = grundgebietstart_x + i * grundgebiet_breite/(n-1)
            y = grundgebietstart_y + j * grundgebiet_breite/(n-1)
            
            #print('i:',i , ',j:', j , ',x:', x ,',y:', y)
            
            _T[i, j] = (T[i+1, j] + T[i-1, j] + T[i, j-1] + T[i, j+1] + 
                        rechte_seite[i,j]
                        ) / 4
            _T[i, j] = T[i, j] + relaxationsfaktor*(_T[i, j] - T[i, j])
    return _T

def einzelschritt_lexikographisch(T,grundgebietstart_x, grundgebietstart_y, grundgebiet_breite, rechte_seite, relaxationsfaktor):

    m, n = T.shape
    _T = np.copy(T)
    # iterate over interior

    for i in range(1, m-1):
        for j in range(1, n-1):
            
            x = grundgebietstart_x + i * grundgebiet_breite/(n-1)
            y = grundgebietstart_y + j * grundgebiet_breite/(n-1)
            
            #print('i:',i , ',j:', j , ',x:', x ,',y:', y)
            
            T[i, j] = (T[i+1, j] + T[i-1, j] + T[i, j-1] + T[i, j+1] + 
                        rechte_seite[i,j]
                        ) / 4
            
            T[i, j] = _T[i, j] + relaxationsfaktor*(T[i, j] - _T[i, j])
    return T

def full_weighting_prolongation(T,homogene_randbedingung):

    m = T.shape[0] * 2 - 1
    n = T.shape[1] * 2 - 1
    
    
    _T = randbedingung_erstellen(n, homogene_randbedingung)
    
    # 4 Situationen werden betrachtet
    # 1. feingitterpunkt = grobgitterpunkt
    for i in range(2, m-1,2):
        for j in range(2, n-1,2):
            
            _T[i,j] = T[int(i/2),int(j/2)]
            
    # 2. feingitterpunkt horizontal zwischen 2 grobgitterpunkten
    for i in range(2, m-2,2):
        for j in range(1, n-1,2):
            
            _T[i,j] = (_T[i,j-1] + _T[i,j+1])/2
            
    # 3. feingitterpunkt vertikal zwischen 2 grobgitterpunkten
    for i in range(1, m-1,2):
        for j in range(2, n-2,2):
            
            _T[i,j] = (_T[i-1,j] + _T[i+1,j])/2
            
    # 4. feingitterpunkt hat 4 diagonale grobgitternachbarn
    for i in range(1, m-1,2):
        for j in range(1, n-1,2):
            
            _T[i,j] = (_T[i-1,j-1] + _T[i+1,j-1] + _T[i-1,j+1] + _T[i+1,j+1])/4
    
    #return np.around(_T,2)
    return _T

def full_weighting_restriktion(T, homogene_randbedingung = 0):
    m = int((T.shape[0] +1) / 2)
    n = int((T.shape[1] +1) / 2)

    _T = randbedingung_erstellen(n, homogene_randbedingung)

    for i in range(1,m-1):
        for j in range(1,n-1):
            _T[i, j] = (T[2*i-1, 2*j-1] + T[2*i-1, 2*j+1] + T[2*i+1, 2*j-1] + T[2*i+1, 2*j+1]
                        + 2*( T[2*i,2*j-1] + T[2*i,2*j+1] + T[2*i-1,2*j] + T[2*i+1, 2*j])
                        + 4* T[2*i,2*j])/16

    return _T

def half_weighting_restriktion(T, homogene_randbedingung = 0):
    m = int((T.shape[0] +1) / 2)
    n = int((T.shape[1] +1) / 2)

    _T = randbedingung_erstellen(n, homogene_randbedingung)

    for i in range(1,m-1):
        for j in range(1,n-1):
            _T[i, j] = (T[2*i,2*j-1] + T[2*i,2*j+1] + T[2*i-1,2*j] + T[2*i+1, 2*j]
                        + 4* T[2*i,2*j])/8

    return _T

def injektion_restriktion(T, homogene_randbedingung = 0):
    m = int((T.shape[0] +1) / 2)
    n = int((T.shape[1] +1) / 2)

    _T = randbedingung_erstellen(n, homogene_randbedingung)

    for i in range(1,m-1):
        for j in range(1,n-1):
            _T[i, j] = T[2*i,2*j]

    return _T

def residuum_berechnen(T,grundgebietstart_x, grundgebietstart_y, grundgebiet_breite, rechte_seite):
    m, n = T.shape

    Residuum = np.zeros((n,n))
    # iterate over interior

    for i in range(1, m-1):
        for j in range(1, n-1):
            
            x = grundgebietstart_x + i * grundgebiet_breite/(n-1)
            y = grundgebietstart_y + j * grundgebiet_breite/(n-1)
            
            #print('i:',i , ',j:', j , ',x:', x ,',y:', y, 'n:', n)
            Residuum[i,j] = rechte_seite[i,j] - ( 
                    4*T[i, j] -T[i+1, j] - T[i-1, j] - T[i, j-1] - T[i, j+1] )
    
    
    return Residuum


def plot(T,grundgebietstart_x, grundgebietstart_y, grundgebiet_breite, exakte_lsg, zykluszahl = -10):
    m, n = T.shape
    
    # x and y axis
    x = np.linspace(grundgebietstart_x, grundgebiet_breite, n)
    y = np.linspace(grundgebietstart_y, grundgebiet_breite, n)
      
    X, Y = np.meshgrid(x, y)
    Z = T
    
    fig, (ax1,ax2,ax3) = plt.subplots(1,3, subplot_kw=dict(projection='3d'))
    
    if(zykluszahl == -10):
        ax1.set_title('Finale Näherungslösung')
    else:
        title = 'Näherungslösung ' + str(zykluszahl)
        ax1.set_title(title)
        
    ax1.plot_surface(X, Y, Z, cmap ='viridis', edgecolor ='green')
    ax1.set_zlim(-6,0)
    ax2.set_title('Exakte Lösung')
    ax2.plot_surface(X, Y, exakte_lsg, cmap ='viridis', edgecolor ='green')
    ax2.set_zlim(-6,0)
    ax3.set_title('Fehler')
    ax3.plot_surface(X, Y, abs(exakte_lsg-Z), cmap ='viridis', edgecolor ='green')
    ax3.set_zlim(0,0.5)
    plt.show()
