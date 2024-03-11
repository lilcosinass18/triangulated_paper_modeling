import matplotlib.pyplot as plt

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
vertices, edges, creases = read_paper_from_file('../txt_info/paper_info2D.txt')

# Визуализация объекта бумаги
fig, ax = plt.subplots()
for edge in edges:
    start, end = edge
    ax.plot([vertices[start][0], vertices[end][0]], [vertices[start][1], vertices[end][1]], 'g-')  # зеленый цвет для ребер

for i, crease in enumerate(creases):
    start, end = crease
    mid_point = ((vertices[start][0] + vertices[end][0]) / 2, (vertices[start][1] + vertices[end][1]) / 2)
    ax.plot([vertices[start][0], vertices[end][0]], [vertices[start][1], vertices[end][1]], 'k--')  # черный пунктир для сгибов
    ax.text(mid_point[0], mid_point[1], str(i), color='purple', ha='center', va='center')

for i, vertex in enumerate(vertices):
    ax.plot(vertex[0], vertex[1], 'bo')  # синий цвет для вершин
    ax.text(vertex[0], vertex[1], str(i), color='red', ha='center', va='bottom')

ax.set_aspect('equal')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Paper Flat Visualization')
plt.grid(True)
plt.savefig('../images/paper_visualization2D.png')
# plt.show()
