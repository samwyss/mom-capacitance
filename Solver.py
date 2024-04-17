from numpy import linspace, meshgrid, ones, zeros, pi, sqrt, floor
from scipy.constants import epsilon_0
from numpy.linalg import solve
import matplotlib.pyplot as plt


class Solver:
    def __init__(self, num_elements_side: int, len_side: float, part: int, phi: float):

        # assign num_elements_side on self
        self.num_elements_side = num_elements_side

        # increment num_elements_side as the following numpy commands operate in terms of vertices
        num_vertices_side = num_elements_side + 1

        # calculate the total number of elements
        num_elements_tot = num_elements_side**2

        # create grid
        # NOTE: for this simple structured grid, all element numbers will be indexed by their upper right node corner
        x = linspace(0, len_side, num_vertices_side)
        y = linspace(0, len_side, num_vertices_side)
        self.X, self.Y = meshgrid(x, y)  # [m]

        # calculate the spacing between vertexes
        self.dist_between_vertex = self.X[0, 0] - self.X[0, 1]  # [m]

        # calculate the area of each element
        self.element_area = self.dist_between_vertex**2  # [m^2]

        # create MoM matrices
        self.b = phi * ones(
            num_elements_tot
        )  # can be assembled without iteration as PHI is constant in simple case
        self.A = zeros((num_elements_tot, num_elements_tot))
        self.c = zeros(num_elements_tot)

        # assemble MoM A matrix based on PART
        if 1 == part:
            for element_m in range(num_elements_tot):
                for element_n in range(num_elements_tot):
                    self.A[element_m, element_n] += self.part_1_a_assembler(
                        element_m, element_n
                    )
        elif 2 == part:
            pass
        elif 3 == part:
            pass
        else:
            raise ValueError("Invalid PART, must be 1, 2 or 3")

    def solve(self):
        # solve linear system of equations
        self.c = solve(self.A, self.b)

        # calculate total charge
        print(self.element_area * sum(self.c) / 1.602e-19)

        plotting_c = (
            self.c.reshape((self.num_elements_side, self.num_elements_side)) / 1.602e-19
        )

        plt.pcolormesh(self.X, self.Y, plotting_c)
        plt.show()

    def part_1_a_assembler(self, element_m: int, element_n: int) -> float:
        if element_m != element_n:
            acc = (
                1
                / (4 * pi * epsilon_0)
                * self.element_area
                / self.calc_element_center_differences(element_m, element_n)
            )
        else:
            acc = 1 / (2 * epsilon_0) * sqrt(self.element_area / pi)
        return acc

    def linear_to_cart_idx(self, element: int) -> tuple[int, int]:
        return int(floor(element / self.num_elements_side)), int(
            element % self.num_elements_side
        )

    def calc_element_center_differences(self, element_m: int, element_n: int) -> float:
        element_m_idx = self.linear_to_cart_idx(element_m)
        element_n_idx = self.linear_to_cart_idx(element_n)

        return sqrt(
            (self.dist_between_vertex * abs(element_m_idx[0] - element_n_idx[0])) ** 2
            + (self.dist_between_vertex * abs(element_m_idx[1] - element_n_idx[1])) ** 2
        )
