import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Функция для чтения информации о бумаге из файла
def read_paper_from_file(filename):
    vertices = []
    edges = []
    creases = []

    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('Vertex:'):
                parts = line.strip().split()
                vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif line.startswith('Edge:'):
                parts = line.strip().split()
                edges.append((int(parts[1]), int(parts[2])))
            elif line.startswith('Crease:'):
                parts = line.strip().split()
                creases.append((int(parts[1]), int(parts[2])))

    return vertices, edges, creases

# Чтение информации о бумаге из файла
vertices, edges, creases = read_paper_from_file('../txt_info/paper_info3D.txt')

# Создание трехмерного графика
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Визуализация ребер
for edge in edges:
    start, end = edge
    x_values = [vertices[start][0], vertices[end][0]]
    y_values = [vertices[start][1], vertices[end][1]]
    z_values = [vertices[start][2], vertices[end][2]]
    ax.plot(x_values, y_values, z_values, 'g-')  # зеленый цвет для ребер

# Визуализация сгибов
for crease in creases:
    start, end = crease
    x_values = [vertices[start][0], vertices[end][0]]
    y_values = [vertices[start][1], vertices[end][1]]
    z_values = [vertices[start][2], vertices[end][2]]
    ax.plot(x_values, y_values, z_values, 'k--')  # черный пунктир для сгибов

# Визуализация вершин
for vertex in vertices:
    ax.scatter(vertex[0], vertex[1], vertex[2], color='blue')  # синий цвет для вершин

# Настройка осей и меток
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Paper Object Visualization')

# Включение поворота
ax.view_init(elev=20, azim=30)

# Отображение графика
plt.show()
