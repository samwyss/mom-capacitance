from scipy import constants as constants

from src.Solver import Solver


def main():
    """
    main driver function
    :return: None
    """

    # tunable constants
    phi = 1  # applied parallel plate voltage
    side_len = 1e-2  # side length of square parallel plate
    num_elements = 10  # number of elements to use per side of square parallel plate
    solver_type = 1  # MoM element formulation number,
    # 1 = circular approximation of square elements
    # 2 = exact square elements
    # 3 = subdomain collocation

    # display empirical results from H. Wintle, “The capacitance of the cube and square plate by random walk
    # methods,” Journal of Electrostatics, vol. 62, no. 51-62, 2004.
    print(
        f"Accepted Empirical Capacitance: {0.36 * 4 * constants.pi * constants.epsilon_0 * side_len} ± {0.01 * 4 * constants.pi * constants.epsilon_0 * side_len} [F]"
    )

    # setup problem and solve
    solver = Solver(num_elements, side_len, solver_type, phi)
    solver.solve()


# call main function after declaration
if __name__ == "__main__":
    main()
