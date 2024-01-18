import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    file_name = 'aaProperties.txt'
    aa_data = pd.read_csv(file_name, sep=';', header=None)

    # Krok 1 - wczytanie danych
    amino_acids = list(aa_data.iloc[:, 0])
    aa_data = aa_data.drop(columns=0)
    aa_data = aa_data.replace(',', '.', regex=True)
    aa_data = aa_data.astype(float)

    # Normalizacja danych
    aa_data = aa_data.apply(lambda x: (x - np.mean(x)) / np.std(x))

    # Krok 2
    cov_matrix = aa_data.cov()

    # Krok 3
    eig_values, eig_vectors = np.linalg.eig(cov_matrix)

    # Sortowanie wektorów własnych
    idx = eig_values.argsort()[::-1]
    eig_values = eig_values[idx]
    eig_vectors = eig_vectors[:, idx]

    # Krok 4
    # Wybór k wektorów własnych
    k = 2
    tmp = eig_vectors[:, :k]
    # Transformacja do nowej przestrzeni (k-wymiarowej)
    aa_data_transformed = aa_data.dot(tmp)

    # Obliczenie wariancji
    variance = eig_values[:2] / eig_values.sum() * 100

    # Rysowanie wykresu
    plt.figure(figsize=(10, 6))
    plt.scatter(aa_data_transformed[0], aa_data_transformed[1])
    plt.xlabel(f'Pierwsza główna składowa ({variance[0]:.2f}% wariancji)')
    plt.ylabel(f'Druga główna składowa ({variance[1]:.2f}% wariancji)')
    for i, aa in enumerate(amino_acids):
        plt.annotate(aa, (aa_data_transformed[0][i], aa_data_transformed[1][i]))
    plt.grid(True)
    plt.show()

    print("Wnioski: ")
    print("Aminokwasy o podobnych właściwościach mają podobne wartości na wykresie.")
    print("Aminokwasy o różnych właściwościach znajdują się w różnych miejscach wykresu, nie są zlokalizowane koło siebie.")
    print("Na przykład aminokwasy: V, L, I, F znajdują się blisko siebie co może sugerować, że mają podobne właściwości.")
