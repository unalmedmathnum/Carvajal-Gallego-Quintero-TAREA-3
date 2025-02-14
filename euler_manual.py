def eulerex(func, ip, iv, fp, h):
    """
    Método de Euler para sistemas de 2 ecuaciones diferenciales ordinarias de primer orden
    Entradas:
        - func: Función con el sistema de EDOs a resolver
        - ip: Punto inicial del intervalo
        - iv: Lista con valor inicial del sistema
        - fp: Punto final del intervalo
        - h: Tamaño de paso
    Salidas:
        - t: Lista conlos puntos donde se aproxima la solución del sistema
        - w0,w1: Listas con las coordenadas del valor de la función en los puntos de t
    """
    w0 = [iv[0]] #Inicialización de los w solución
    w1 = [iv[1]]
    t = [ip] #Diccionario con los puntos donde se aproxima la solución
    
    p = ip
    y = iv.copy() 
    
    while p < fp:
        
        dy = func(p, y)
        y = [y[i] + h * dy[i] for i in range(len(y))] #Creación del siguiente w
        p += h #Actualización de paso

        t.append(p) #Guardar resultados
        w0.append(y[0])
        w1.append(y[1])
        
    return t, w0, w1