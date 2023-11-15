import matplotlib.pyplot as plt
import numpy as np

def cubic_spline_natural(x, y):
    """
    Calcula os coeficientes dos polinômios cúbicos naturais para interpolação.

    Args:
        x (list): Lista de coordenadas x dos pontos de interpolação.
        y (list): Lista de coordenadas y dos pontos de interpolação.

    Returns:
        tuple: Coeficientes dos polinômios cúbicos naturais (a, b, c, d).
    """
    n = len(x)
    h = {k: x[k+1] - x[k] for k in range(n - 1)}

    A = np.zeros((n, n))
    for i in range(1, n - 1):
        A[i, i-1] = h[i-1]
        A[i, i] = 2 * (h[i-1] + h[i])
        A[i, i+1] = h[i]

    A[0, 0] = 1
    A[-1, -1] = 1

    B = np.zeros(n)
    for k in range(1, n - 1):
        B[k] = 3 * ((y[k+1] - y[k]) / h[k] - (y[k] - y[k-1]) / h[k-1])

    c = np.linalg.solve(A, B)

    a = y
    b = np.zeros(n-1)
    d = np.zeros(n-1)
    for k in range(n-1):
        b[k] = (1/h[k]) * (a[k+1] - a[k]) - (h[k]/3) * (2*c[k] + c[k+1])
        d[k] = (c[k+1] - c[k]) / (3 * h[k])

    return a, b, c, d

def evaluate_spline(x, a, b, c, d, xi):
    """
    Avalia o valor da spline cúbica no ponto xi.

    Args:
        x (list): Lista de coordenadas x dos pontos de interpolação.
        a (list): Coeficientes a dos polinômios cúbicos naturais.
        b (list): Coeficientes b dos polinômios cúbicos naturais.
        c (list): Coeficientes c dos polinômios cúbicos naturais.
        d (list): Coeficientes d dos polinômios cúbicos naturais.
        xi (float): Ponto onde a spline cúbica é avaliada.

    Returns:
        float: Valor da spline cúbica no ponto xi.
    """
    k = 0
    while x[k+1] < xi:
        k += 1

    dx = xi - x[k]
    result = a[k] + b[k] * dx + c[k] * dx**2 + d[k] * dx**3
    return result

# Entrada dos valores via teclado
n = int(input("Digite a quantidade de pontos de interpolação: "))
x_values = list(map(float, input("Digite os valores de x separados por espaço: ").split()))
y_values = list(map(float, input("Digite os valores de y separados por espaço: ").split()))
z_value = float(input("Digite o valor de z: "))
m = int(input("Digite a quantidade de pontos m para calcular as imagens das splines em D: "))

# Verifica se o número de pontos de interpolação é igual
if len(x_values) != n or len(y_values) != n:
    print("Erro: O número de pontos de interpolação não é igual.")
else:
    a, b, c, d = cubic_spline_natural(x_values, y_values)

    # Calcular si(z)
    si_z = evaluate_spline(x_values, a, b, c, d, z_value)
    print(f"Para z = {z_value}, si(z) = {si_z}")

    # Imprimir valores xi e suas respectivas imagens si(x) para o conjunto de m pontos em D
    print("Valores de xi e suas respectivas imagens si(x) para o conjunto de m pontos em D:")
    for xi in np.linspace(x_values[0], x_values[-1], m):
        si_xi = evaluate_spline(x_values, a, b, c, d, xi)
        print(f"xi = {xi}, si(xi) = {si_xi}")

    # Gerar pontos igualmente espaçados em D para plotagem
    plot_points = np.linspace(x_values[0], x_values[-1], m)

    # Calcular as imagens das splines em D
    spline_images = [evaluate_spline(x_values, a, b, c, d, xi) for xi in plot_points]

    # Plotar as splines em D
    plt.scatter(x_values, y_values, color='red', label='Pontos de Interpolação')
    plt.plot(plot_points, spline_images, label='Splines em D')
    plt.legend()
    plt.savefig('spline.png')
