from Solver import Solver

PHI = 1
SIDE_LEN = 1e-2


def main():
    # part 1
    print("Part 1")
    solver_1_10 = Solver(10, SIDE_LEN, 1, PHI)
    solver_1_20 = Solver(20, SIDE_LEN, 1, PHI)
    solver_1_30 = Solver(30, SIDE_LEN, 1, PHI)
    solver_1_10.solve()
    solver_1_20.solve()
    solver_1_30.solve()
    print("\n")

    # part 2
    print("Part 2")
    solver_2_10 = Solver(10, SIDE_LEN, 2, PHI)
    solver_2_20 = Solver(20, SIDE_LEN, 2, PHI)
    solver_2_30 = Solver(30, SIDE_LEN, 2, PHI)
    solver_2_10.solve()
    solver_2_20.solve()
    solver_2_30.solve()
    print("\n")


# call main function after declaration
if __name__ == "__main__":
    main()
