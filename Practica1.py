#Práctica 1
#Alumno: WALID SABHI

from Bio import SeqIO
import numpy as np
from Bio.Align import Alignment,substitution_matrices
from Bio.Seq import Seq

"""
Paso 1: importar las secuencias de los archivos .fasta
"""
#importar las dos secuencias de los archivos .fasta
iterador_reg=SeqIO.parse(r"https://github.com/flooki10/Alineamiento_Pareado/blob/main/seq1.fasta","fasta")
iterador_reg2=SeqIO.parse(r"https://github.com/flooki10/Alineamiento_Pareado/blob/main/seq2.fasta","fasta")
#imprimimos las dos secuencias 
secuencia_1=next(iterador_reg)
secuencia_2=next(iterador_reg2)
print(secuencia_1)
print(secuencia_2)



"""
Paso 2: creamos la función needleman wunsch
"""
#Ahora vamos a hacer el alineamiento pareado utilizando las dos secuencias
#NeedlemanWunsh para el alineamiento global
#tenemos que calcular las matrices N y T, empleando un valor de penalización por gap de -8
def needleman_wunsch(seq_query, seq_target, score_matrix,  gap_pen):

    #creamos las matices N y de T de ceros  
    len_query = len(seq_query)
    len_target = len(seq_target)
    N = np.zeros((len_query + 1, len_target + 1), dtype=int)
    T = np.zeros((len_query + 1, len_target + 1), dtype=int)

    #ahora ya como tenemos las matrices creadas y ya hemos rellenado la fila y columna estandar pasamos a rellenar nuestras matices
    for i in range(1,len(seq_query)+1):
        for j in range(1,len(seq_target)+1):
            #aqui lo que he hecho es recorrer la primera columna y la primera fila de las matrices de N y T y colocar los valores estandares de dirección y de gap
            N[i][0]=gap_pen*i
            T[i][0]=2

            N[0][j]=gap_pen*j
            T[0][j]=3

            #trabajamos primero sobre la matriz N
            diagonal= N[i-1][j-1]+score_matrix[seq_query[i-1]][seq_target[j-1]]
            arriba= N[i-1][j]+gap_pen
            izquierda= N[i][j-1]+gap_pen

            valor_max= max(diagonal,arriba,izquierda)

            N[i][j]=valor_max
            if valor_max == diagonal:
                T[i][j]= 1
            elif valor_max == arriba:
                T[i][j]=2
            else:
                T[i][j]=3
    return N,T



gap_pen=-8
#Metodo 1
#definimos la matriz de sustitución
score_matrix={
    'A': {'A': 9, 'T': -2, 'G': -3, 'C': -4},
    'T': {'A': -2, 'T': 7, 'G': -5, 'C': -2},
    'G': {'A': -3, 'T': -5, 'G': 10, 'C': -1},
    'C': {'A': -4, 'T': -2, 'G': -1, 'C': 7}
                }

#Metodo 2
#puedo también que probar hacerla con Arrays (TEMA2_B)

#imprimimos la matiz de deedlman
N,T= needleman_wunsch(secuencia_1, secuencia_2, score_matrix,  gap_pen)
print(N)
print("\n")
print(T)

"""
#PASO3: Escritura del alineamiento

Metodo 2
#PASO3: Escritura del alineamiento
alineador=Align.PairwiseAligner(mode='global')
alineamientos=alineador.align(secuencia_2,secuencia_1)
alineamiento=alineamientos[0]
print(alineamiento.shape)
print(alineamiento)

#sacamos la puntuación obtenida en la matriz N obtenida
puntuacion= N[len(secuencia_1)][len(secuencia_2)]
print("Puntuación del alineamiento: ",puntuacion)

"""

#vamos a generar el alineamiento a partir de las matrices N y T
def generar_aliniamiento(seq_query, seq_target, N, T):
    alignment_query = ""
    alignment_target = ""
    i, j = len(seq_query), len(seq_target)
    
    while i > 0 or j > 0:
        if T[i][j] == 1:  #si es match o mismatch
            alignment_query += seq_query[i-1] 
            alignment_target += seq_target[j-1] 
            i -= 1
            j -= 1
        #valor en la matriz T es 2, es un gap en la secuencia target
        elif T[i][j] == 2:
            alignment_query = seq_query[i - 1] + alignment_query
            alignment_target = "-" + alignment_target
            i -= 1
        #valor en la matriz T es 3, es un gap en la secuencia query
        else:
            alignment_query = "-" + alignment_query
            alignment_target = seq_target[j - 1] + alignment_target
            j -= 1

    return alignment_query, alignment_target


alignment_query,alignment_target=generar_aliniamiento(secuencia_1, secuencia_2, N, T)
print(alignment_query,alignment_target)


def Alineamiento(secuencias):
    #quitamos - de las secuencias
    sequences = [line.replace("-", "") for line in secuencias]
    #print(sequences)

    coordinates = Alignment.infer_coordinates(secuencias)
    #print(coordinates)

    alignment = Alignment(sequences, coordinates)
    print(alignment)

    #obtención del puntuacion del alineamiento que también lo podemos sacar de 2 maneras
    #una seria obtener el mayor numero de la matriz usando len el secundo metodo seria uasr max de la matriz
    puntuacion= N[len(secuencia_1)][len(secuencia_2)]
    print("Puntuación del alineamiento: ",puntuacion)



#secuencias tiene que ser una lista de las dos secuencias 
secuencias  =[alignment_query,alignment_target]
Alineamiento(secuencias)

print("\n")

"""
Paso 4: Uso de pipline para la comparación y análisis de similitud entre dos secuencias,
        una de ADN y otra de proteina obtenidas de pubmed

        puntos a discutir:
        • El origen y el formato de las secuencias.
        • La naturaleza biológica de las secuencias (especie, gen, función, etc). Razona por qué lo has elegido.
        • ¿Qué partes se han conservado más evolutivamente y cuáles menos? ¿A qué crees que se deben las
        diferencias a nivel de conservación entre diferentes zonas?
        • Utiliza Biopython para transformar las secuencias de ADN que has seleccionado a secuencias de aminoácidos
        (proteínas). Ejecuta el pipeline para realizar su alineamiento después de la traducción. Explica cómo has
        realizado el proceso de traducción e incluye y analiza el resultado obtenido tras el alineamiento.


        #usando NIH y uniprot
        #gen BRCA1 en Homosapiens(humano) y Musculus(ratones), asociado con el cáncer de mama y ovario
        #vamos a elegir una pareja de secuencias de ADN en Homosapiens y musculus y otra de proteinas del gen BRCA1
        #pareja de ADN(NIH): https://www.ncbi.nlm.nih.gov/nuccore/OR780020. / https://www.ncbi.nlm.nih.gov/nuccore/EU349657.1
        #pareja de proteina(Uniprot): https://www.uniprot.org/uniprotkb/P38398/entry  / https://www.uniprot.org/uniprotkb/P48754/entry
        hemos descargados las secuencias en archivos .fasta para facilitar su uso luego

"""


#importar las dos secuencias de los archivos .fasta
iterador_reg1=SeqIO.parse(r"C:\Users\walid\OneDrive\Desktop\bioinformatica\P1_SABHI_WALID\seq_ADN_BRCA1_humano.fasta","fasta")
iterador_reg2=SeqIO.parse(r"C:\Users\walid\OneDrive\Desktop\bioinformatica\P1_SABHI_WALID\seq_ADN_BRCA1_raton.fasta","fasta")
iterador_reg3=SeqIO.parse(r"C:\Users\walid\OneDrive\Desktop\bioinformatica\P1_SABHI_WALID\seq_PROTEINA_BRCA1_humano.fasta","fasta")
iterador_reg4=SeqIO.parse(r"C:\Users\walid\OneDrive\Desktop\bioinformatica\P1_SABHI_WALID\seq_PROTEINA_BRCA1_raton.fasta","fasta")
#imprimimos las dos secuencias 
seq_ADN_humano=next(iterador_reg1)
seq_ADN_raton=next(iterador_reg2)
seq_PROTEINA_humano=next(iterador_reg3)
seq_PROTEINA_raton=next(iterador_reg4)

print(seq_ADN_humano)
print(seq_ADN_raton)
print(seq_PROTEINA_humano)
print(seq_PROTEINA_raton)


print("\n")

#como que en este caso tenemos una secuencia de proteina debemos cargar la matriz de sustitucion de proteinas
pam250 = substitution_matrices.load("PAM250")
#print(pam250)

#
def needleman_wunsch_2(seq_query, seq_target, score_matrix, pam250,  gap_pen):

    #creamos las matices N y de T de ceros  
    len_query = len(seq_query)
    len_target = len(seq_target)
    N = np.zeros((len_query + 1, len_target + 1), dtype=int)
    T = np.zeros((len_query + 1, len_target + 1), dtype=int)

    #aqui lo que he hecho es recorrer la primera columna y la primera fila de las matrices de N y T y colocar los valores estandares de dirección y de gap
    for i in range(1,len(seq_query)+1):
        N[i][0]=gap_pen*i
        T[i][0]=2

    for j in range(1,len(seq_target)+1):
        N[0][j]=gap_pen*j
        T[0][j]=3

    #ahora ya como tenemos las matrices creadas y ya hemos rellenado la fila y columna estandar pasamos a rellenar nuestras matices
    for i in range(1,len(seq_query)+1):
        for j in range(1,len(seq_target)+1):
            if (seq_query[i - 1] in score_matrix) and (seq_target[j - 1] in score_matrix):
                score = score_matrix
            else:
                score = pam250
                
            diagonal = N[i - 1][j - 1] + score[seq_query[i - 1]][seq_target[j - 1]]
            arriba= N[i-1][j]+gap_pen
            izquierda= N[i][j-1]+gap_pen

            valor_max= max(diagonal,arriba,izquierda)

            N[i][j]=valor_max
            if valor_max == diagonal:
                T[i][j]= 1
            elif valor_max == arriba:
                T[i][j]=2
            else:
                T[i][j]=3
    return N,T

"""
comparación de la pareja de secuencias de ADN
"""

print("Comparación de la pareja de secuencias de ADN")
N,T=needleman_wunsch(seq_ADN_humano, seq_ADN_raton, score_matrix,  gap_pen)
#imprimir el resultado de la funcion
print("Matriz N","\n",N)
print("\n")
print("Matriz T","\n",T)

#alineamiento
x,y=generar_aliniamiento(seq_ADN_humano, seq_ADN_raton, N, T)
print(x,y)

Alineamiento([x,y])

"""
comparación de la pareja de secuencias de proteinas 
"""
print("\n")
print("Comparación de lapareja de secuencias de proteinas")
N,T=needleman_wunsch_2(seq_PROTEINA_humano,seq_PROTEINA_raton,score_matrix,pam250,gap_pen)
#imprimir el resultado de la funcion
print("Matriz N","\n",N)
print("\n")
print("Matriz T","\n",T)

#alineamiento
x,y=generar_aliniamiento(seq_PROTEINA_humano, seq_PROTEINA_raton, N, T)
print(x,y)

Alineamiento([x,y])


""" Pregunta !!!
#local
alineador_local=Align.PairwiseAligner(mode='local')
alineamientos_locales=alineador_local.align(secuencia_2,secuencia_1)
alineamiento_local=alineamientos_locales[0]
print(alineamiento_local.shape)
print(alineamiento_local)
"""

