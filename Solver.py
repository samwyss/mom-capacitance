from numpy import linspace, meshgrid, ones, zeros, pi, sqrt, floor, log
from scipy.constants import epsilon_0, e
from numpy.linalg import solve
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


class Solver:
    def __init__(self, num_elements_side: int, len_side: float, part: int, phi: float):
        """
        Solver constructor
        :param num_elements_side: number of elements per side of grid
        :param len_side: length of the side in meters
        :param part: method by which to calculate A_mn from textbook
        :param phi: potential of the plate
        """

        # instantiate self.capacitance and self.total_charge to zero
        self.capacitance = 0.0
        self.total_charge = 0.0

        # assign num_elements_side on self
        self.num_elements_side = num_elements_side

        # assign phi on self
        self.phi = phi

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
        self.dist_between_vertex = self.X[0, 1] - self.X[0, 0]  # [m]

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
            for element_m in range(num_elements_tot):
                for element_n in range(num_elements_tot):
                    self.A[element_m, element_n] += self.part_2_a_assembler(
                        element_m, element_n
                    )
        elif 3 == part:
            for element_m in range(num_elements_tot):
                for element_n in range(num_elements_tot):
                    pass
        else:
            raise ValueError("Invalid PART, must be 1, 2 or 3")

    def solve(self):
        """
        solves the system of equations set up in the Solver constructor
        :return: None
        """
        # solve linear system of equations
        self.c = solve(self.A, self.b)

        # calculate total charge
        self.total_charge = self.element_area * sum(self.c)

        # print total charge out
        print(
            f"{self.num_elements_side}x{self.num_elements_side} Total Charge: {self.total_charge} [C]"
        )

        # calculate capacitance
        self.capacitance = self.total_charge / self.phi

        # print out capacitance
        print(
            f"{self.num_elements_side}x{self.num_elements_side} Capacitance: {self.capacitance} [F]"
        )

        # reshape results for plotting
        reshaped_charge = self.c.reshape(
            (self.num_elements_side, self.num_elements_side)
        )

        # plot charge distribution
        plt.rcParams["text.usetex"] = True
        fig, ax = plt.subplots(dpi=300, figsize=(5, 4))
        mesh = ax.pcolormesh(
            self.X,
            self.Y,
            reshaped_charge,
            cmap="jet",
            norm=LogNorm(vmin=1e-9, vmax=2e-8),
        )
        cbar = plt.colorbar(mesh)
        fig.show()

    def part_1_a_assembler(self, element_m: int, element_n: int) -> float:
        """
        calculates A_mn for MoM matrices for part 1
        :param element_m: element m
        :param element_n: element n
        :return: A_mn
        """
        if element_m != element_n:
            amn = (
                1
                / (4 * pi * epsilon_0)
                * self.element_area
                / self.calc_element_center_differences(element_m, element_n)
            )
        else:
            amn = 1 / (2 * epsilon_0) * sqrt(self.element_area / pi)
        return amn

    def linear_to_cart_idx(self, element: int) -> tuple[int, int]:
        """
        converts a linear index to a cartesian grid
        :param element: linearly indexed element
        :return: cartesian coordinates of the linear index
        """
        return int(element % self.num_elements_side), int(
            floor(element / self.num_elements_side)
        )

    def calc_element_center_differences(self, element_m: int, element_n: int) -> float:
        """
        calculates the difference in center point location of element m and element n which are represented by their
        top left nodes
        :param element_m: element m
        :param element_n: element n
        :return: distance between the center points of elements m and n
        """
        element_m_idx = self.linear_to_cart_idx(element_m)
        element_n_idx = self.linear_to_cart_idx(element_n)

        return sqrt(
            (self.dist_between_vertex * abs(element_m_idx[0] - element_n_idx[0])) ** 2
            + (self.dist_between_vertex * abs(element_m_idx[1] - element_n_idx[1])) ** 2
        )

    def part_2_a_assembler(self, element_m: int, element_n: int) -> float:
        """
        calculates A_mn for MoM matrices for part 2
        :param element_m: element m
        :param element_n: element n
        :return: A_mn
        """

        # ensure accumulator is initialized to zero
        amn = 0

        # get centerpoint of elements m and n
        xm, ym = self.get_element_center_point(element_m)
        xn, yn = self.get_element_center_point(element_n)

        # get local subgrid of "prime" variable
        xp = [
            xn - self.dist_between_vertex * 0.5,
            xn + self.dist_between_vertex * 0.5,
        ]
        yp = [
            yn - self.dist_between_vertex * 0.5,
            yn + self.dist_between_vertex * 0.5,
        ]

        # accumulate on amn
        for i in range(2):
            for j in range(2):
                # determine the sign of the term
                sign = 1 if (i + j) % 2 == 0 else -1

                # calculate R
                r = sqrt((xm - xp[i]) ** 2 + (ym - yp[j]) ** 2)

                # initialize values in case the log term evaluates to zero
                term_1 = 0
                term_2 = 0

                # calculate first term contributions
                if ((ym - yp[j]) + r) != 0:
                    term_1 = (xm - xp[i]) * log((ym - yp[j]) + r)

                # calculate second term contributions
                if ((xm - xp[i]) + r) != 0:
                    term_2 = (ym - yp[j]) * log((xm - xp[i]) + r)

                # accumulate on amn
                amn += sign * (term_1 + term_2)

        # multiply by constant and return
        return amn * 1 / (4 * pi * epsilon_0)

    def get_element_center_point(self, element: int):
        """
        calculates the center point location of a linearly indexed element
        :param element: linearly indexed element
        :return: center point location of linearly indexed element
        """
        coords = self.linear_to_cart_idx(element)

        return (
            self.X[coords[0], coords[1]] + self.dist_between_vertex * 0.5,
            self.Y[coords[0], coords[1]] + self.dist_between_vertex * 0.5,
        )
